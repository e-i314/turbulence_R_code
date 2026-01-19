
library(archive)
library(httr2)
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(data.table)
library(signal)
library(purrr)
library(stringr)
library(tibble)
library(tidyr)

# Download the ZIP from Zenodo ---
zenodo_record_id <- "17988010"  
zip_name <- "turbulence_dataset.zip"

zip_url <- sprintf(
  "https://zenodo.org/records/%s/files/%s?download=1",
  zenodo_record_id, zip_name
)

dest_zip <- file.path(tempdir(), zip_name)

options(timeout = 1200)  # 20min
download.file(zip_url, destfile = dest_zip, mode = "wb", method = "libcurl")

# extract to  C:\Users/yourfolder/document/zenodo_turbulence_dataset/

extract_dir <- file.path(path.expand("~"), "zenodo_turbulence_dataset")
dir.create(extract_dir, showWarnings = FALSE)
unzip(dest_zip, exdir = extract_dir, overwrite = TRUE)

dataset_root <- extract_dir

# --- ULSA PRO (1 Hz) ---
pole_1hz_path <- file.path(dataset_root, "raw", "pole_ULSA_PRO",
                           "ULSA_PRO_1hz_20250114_20251207.csv")
pole_1hz <- readr::read_csv(pole_1hz_path, show_col_types = FALSE)



# ----------------------------
# 0. to_1Hz
# ----------------------------




setDT(pole_1hz)  # df を data.table化（コピーしない）
pole_1hz[, time := as.POSIXct(`_time`)]
pole_1hz[, time_1s := as.POSIXct(floor(as.numeric(time)), origin="1970-01-01", tz=attr(time,"tzone"))]

# 循環平均（deg）
circ_mean_deg_dt <- function(theta){
  th <- theta * pi/180
  (atan2(mean(sin(th), na.rm=TRUE), mean(cos(th), na.rm=TRUE)) * 180/pi + 360) %% 360
}

df <- pole_1hz[
  ,
  .(
    windspd      = median(windspd, na.rm=TRUE),
    windspd_head = median(windspd_head, na.rm=TRUE),
    temp = median(temp, na.rm = TRUE),
    winddir      = circ_mean_deg_dt(winddir)
  ),
  by = time_1s
][order(time_1s)]

setnames(df, "time_1s", "_time")


# ----------------------------
# 1. MagHeading
# ----------------------------

#
df <- df %>%
  mutate(
    winddir_mag = (winddir - 46 + 360) %% 360
  )


# ----------------------------
# 2. time & 10min window
# ----------------------------
df <- df %>%
  mutate(
    time = as.POSIXct(`_time`, tz = "UTC"),
    block_10min = floor_date(time, "10 minutes")
  )

# ----------------------------
# 3. wind direction average each 10min windows
# ----------------------------
mean_dir <- df %>%
  group_by(block_10min) %>%
  summarise(
    mean_dir_rad = atan2(
      mean(sin(winddir * pi / 180), na.rm = TRUE),
      mean(cos(winddir * pi / 180), na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    mean_dir = (mean_dir_rad * 180 / pi + 360) %% 360
  )


# ----------------------------
# 5. Projection to wind direction
# ----------------------------
# 風向差（ラジアン）
df2 <- df %>%
  left_join(mean_dir, by = "block_10min") %>%
  mutate(
    # 各時刻の風向と10分平均風向の差（deg→rad）
    dtheta = (winddir - mean_dir) * pi/180,
    
    # 10分平均風向方向への射影（割り算しない）
    spd_proj = windspd * cos(dtheta)
  ) %>%
  arrange(time) %>%
  group_by(block_10min) %>%
  mutate(
    d1 = spd_proj - lag(spd_proj, 1),
    d2 = spd_proj - lag(spd_proj, 2),
    d5 = spd_proj - lag(spd_proj, 5)
  ) %>%
  ungroup()

# ----------------------------
# 6. du（1s, 2s, 5s）
#    
# ----------------------------
df2 <- df2 %>%
  arrange(time) %>%
  group_by(block_10min) %>%
  mutate(
    dU_1s = spd_proj - lag(spd_proj, 1),
    dU_2s = spd_proj - lag(spd_proj, 2),
    dU_5s = spd_proj - lag(spd_proj, 5),
    dU_head =windspd_head - lag(windspd_head, 1)
  ) %>%
  ungroup()


# ----------------------------
# 5. ΔU* = (ΔU - median) / b_k
#    b_k = mean(|ΔU|)
# ----------------------------
df_norm <- df2 %>%
  group_by(block_10min) %>%
  mutate(
    center_1s = median(dU_1s, na.rm = TRUE),
    b_k_1s    = mean(abs(dU_1s - center_1s), na.rm = TRUE),
    dU_std = (dU_1s - center_1s) / b_k_1s,
    
    center_head = median(dU_head, na.rm = TRUE),
    b_head    = mean(abs(dU_head - center_head), na.rm = TRUE),
    dU_std_head = (dU_head - center_head) / b_head    
    
    
  ) %>%
  ungroup() %>%
  dplyr::filter(is.finite(dU_std))

q <- quantile(df_norm$b_head, probs = c(0.2, 0.8), na.rm = TRUE)

df_norm <- df_norm %>%
  mutate(
    b_group = case_when(
      b_head <= q[1] ~ "Low b",
      b_head >= q[2] ~ "High b",
      TRUE        ~ "Mid b"
    )
  )


# ----------------------------
# 6. PDF
# ----------------------------



# 1) 10分窓ごとの平均風速を計算
win_keep <- df_norm %>%
  group_by(block_10min) %>%
  summarise(mean_windspd = mean(windspd, na.rm = TRUE), .groups = "drop") %>%
#  dplyr::filter(mean_windspd >= 1) %>%
  select(block_10min)

# 2) 図1に使うデータ（= 条件を満たす窓の全サンプル）を抽出
df_fig1 <- df_norm %>%
  inner_join(win_keep, by = "block_10min") %>%
  dplyr::filter(is.finite(dU_std), is.finite(dU_std_head))

# （任意）窓数と点数の確認
n_win <- n_distinct(df_fig1$block_10min)
n_pts <- nrow(df_fig1)
cat("n_win =", n_win, " / n_pts =", n_pts, "\n")

# Draw fig_1a
ggplot(df_fig1) +
  geom_density(
    aes(dU_std, colour = "Along-wind"),
    size = 1,
    n = 4096,
    bw = 0.3
  ) +
  geom_density(
    aes(dU_std_head, colour = "Head-axis"),
    size = 1,
    n = 4096,
    bw = 0.3
  ) +
  scale_colour_manual(
    values = c(
      "Along-wind" = "black",
      "Head-axis"  = "#0072B2"
    )
  ) +
  scale_y_log10() +
  labs(
    x = expression(Delta*U^"*"),
    y = "Probability density",
    colour = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    axis.title = element_text(size = 9 * 1.3),
    axis.text  = element_text(size = 9 * 1.3),
    legend.text = element_text(size = 9 * 1.3)
  ) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(1e-4, 1))

# Draw fig_1b
ggplot(df_norm, aes(dU_head, colour = b_group)) +
  geom_density(size = 1) +
  scale_y_log10() +
  scale_colour_discrete(name = NULL)+
  labs(
    x = expression(Delta*U[head] (m/s)),
    y = "Probability density"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 9 * 1.3),
    axis.text  = element_text(size = 9 * 1.3),
    legend.text = element_text(size = 9 * 1.3)
  )+
  xlim(-10,10)+
  coord_cartesian(xlim = c(-5, 5), ylim = c(1e-3, 1e3))




df2_std <- df2%>%
  mutate(
    time = as.POSIXct(time),
    windspd = as.numeric(windspd),
    winddir = as.numeric(winddir),
    u = - windspd * sin(winddir * pi / 180),
    v = - windspd * cos(winddir * pi / 180)
  ) %>%
  arrange(time) %>%
  mutate(window = floor_date(time, "10 minutes"))




project_velocity <- function(u, v, theta_deg){
  th <- theta_deg * pi / 180
  u * cos(th) + v * sin(th)
}

laplace_mle_b <- function(x, min_n = 50){
  x <- x[is.finite(x)]
  if(length(x) < min_n) return(NA_real_)
  mu <- median(x)
  mean(abs(x - mu))
}

angles <- seq(0, 350, by = 10)
min_n  <- 50

b_list <- list()   

# Take 60min for caliculate
for(w in unique(df2_std$window)){
  
  dat <- df2_std %>%
    dplyr::filter(window == w) %>%
    arrange(time)
  
  # 点数が足りない窓はスキップ
  if(sum(is.finite(dat$u) & is.finite(dat$v)) < (min_n + 1)){
    next
  }
  
  for(a in angles){
    
    uproj <- project_velocity(dat$u, dat$v, a)
    du    <- diff(uproj)
    
    b_val <- laplace_mle_b(du, min_n = min_n)
    
    b_list[[length(b_list) + 1]] <- tibble(
      window = w,
      angle  = a,
      b      = b_val
    )
  }
}

b_df <- bind_rows(b_list)

# 「各10分窓内」で角度を振ったときのばらつき

win_summary <- b_df %>%
  group_by(window) %>%
  summarise(
    n_angle = sum(is.finite(b)),
    med_b   = median(b, na.rm = TRUE),
    mean_b  = mean(b, na.rm = TRUE),
    sd_b    = sd(b, na.rm = TRUE),
    cv_b    = sd_b / mean_b,                     # 角度方向の変動係数
    b_min   = min(b, na.rm = TRUE),
    b_max   = max(b, na.rm = TRUE),
    ratio_maxmin = b_max / b_min,                # 角度方向 max/min（窓内）
    rel_range = (b_max - b_min) / med_b,         # 窓内中央値に対する相対レンジ
    .groups = "drop"
  ) %>%
  dplyr::filter(is.finite(cv_b), is.finite(ratio_maxmin), n_angle >= 10)

# 要約（中央値と95%点など）
summ_stats <- win_summary %>%
  summarise(
    n_win = n(),
    cv_med  = median(cv_b, na.rm = TRUE),
    cv_p95  = quantile(cv_b, 0.95, na.rm = TRUE),
    ratio_med = median(ratio_maxmin, na.rm = TRUE),
    ratio_p95 = quantile(ratio_maxmin, 0.95, na.rm = TRUE),
    relrng_med = median(rel_range, na.rm = TRUE),
    relrng_p95 = quantile(rel_range, 0.95, na.rm = TRUE)
  )

summ_stats

# Draw fig_2
ggplot(b_df, aes(x = factor(angle), y = b)) +
  geom_boxplot() +
  labs(
    x = "Projection angle from nose direction (deg)",
    y = "b (m/s)"
  ) +
  theme_bw()+
  theme(    axis.title = element_text(size = 9 * 1.3),
            axis.text  = element_text(size = 9 * 1.3),
            axis.text.x = element_text(
              angle = 90,
              vjust = 0.5,
              hjust = 1) )+
  scale_y_log10()



q95 <- quantile(b_df$b, 0.95, na.rm = TRUE)

t_center_unix <- b_df %>%
  dplyr::filter(is.finite(b), b >= q95) %>%
  dplyr::arrange(b) %>%
  dplyr::slice(1) %>%
  dplyr::pull(window)
t_start <- as.POSIXct(t_center_unix, origin = "1970-01-01", tz = "Asia/Tokyo")
t_end    <- t_start + minutes(10)

# ---------- parameter ----------
time_col <- "time"  
b_col    <- "b"     
x_col    <- "windspd_head"     # 分布比較したい系列（例：u 成分）
Delta_t  <- 1     # Δt
win      <- minutes(10)
df10 <- df2 %>% dplyr::filter(time >= t_start,
                          time <= t_end)

# ---------- 固定Δt の増分 δx(t;Δt) ----------
dt0 <- median(as.numeric(diff(df10$time), units = "secs"), na.rm = TRUE) # 実質サンプリング間隔
k   <- max(1L, as.integer(round(Delta_t / dt0)))                        # Δt に相当するラグ点数

x  <- df10[[x_col]]
dx <- x - dplyr::lag(x, k)
dx <- dx[!is.na(dx)]

# ---------- 当てはめ（正規：平均・標準偏差，ラプラス：median と MAD(平均絶対偏差)） ----------
muN <- mean(dx); sdN <- sd(dx)
muL <- median(dx)
bL  <- mean(abs(dx - muL))  # ラプラス尺度（ユーザ変数 b と混同しないため bL とする）

qlaplace <- function(p, location = 0, scale = 1){
  location + scale * ifelse(p < 0.5, log(2*p), -log(2*(1-p)))
}

# ---------- Q-Q  ----------
n  <- length(dx)
pp <- ppoints(n)
samp <- sort(dx)

qq_norm <- tibble(theo = qnorm(pp, muN, sdN), samp = samp)
qq_lap  <- tibble(theo = qlaplace(pp, muL, bL), samp = samp)

qqline_par <- function(theo, samp){
  qy <- quantile(samp, c(0.25, 0.75))
  qx <- quantile(theo, c(0.25, 0.75))
  a  <- diff(qy) / diff(qx)
  b  <- qy[1] - a * qx[1]
  list(slope = a, intercept = b)
}
ln1 <- qqline_par(qq_norm$theo, qq_norm$samp)
ln2 <- qqline_par(qq_lap$theo,  qq_lap$samp)




# dx
n  <- length(dx)
pp <- ppoints(n)
samp <- sort(dx)

# normal
muN <- mean(dx)
sdN <- sd(dx)
qq_norm <- tibble(
  theo = qnorm(pp, muN, sdN),
  samp = samp
)

# laplace
muL <- median(dx)
bL  <- mean(abs(dx - muL))
qlaplace <- function(p, mu, b){
  mu + b * ifelse(p < 0.5, log(2*p), -log(2*(1-p)))
}
qq_lap <- tibble(
  theo = qlaplace(pp, muL, bL),
  samp = samp
)


# Caliculate fig_3
df_ts <- df10 %>%
  mutate(
    time = with_tz(time, "Asia/Tokyo")
  ) %>%
  select(time, windspd, windspd_head) %>%
  pivot_longer(
    cols = c(windspd, windspd_head),
    names_to = "variable",
    values_to = "value"
  )


# Draw fig_3a
ggplot(df_ts, aes(x = time, y = value)) +
  geom_line(linewidth = 0.3, colour = "black") +
  facet_wrap(
    ~ variable,
    ncol = 1,
    scales = "free_y",
    labeller = as_labeller(c(
      windspd      = "Wind speed",
      windspd_head = "Wind speed head"
    ))
  ) +
  labs(
    x = "Time",
    y = "Wind speed (m/s)"
  ) +
  theme_bw(base_size = 9) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 9 * 1.3),  # ← facet文字を1.3倍
    panel.spacing = unit(2, "mm"),
    axis.title = element_text(size = 9 * 1.3),
    axis.text  = element_text(size = 9 * 1.3)
  )


#Draw fig_3b
ggplot() +
  # empirical histogram (very light)
  geom_histogram(
    aes(dx, after_stat(density)),
    bins = 60,
    fill = "grey80",
    colour = NA
  ) +
  # Normal PDF
  stat_function(
    aes(colour = "Normal"),
    fun = dnorm,
    args = list(mean = muN, sd = sdN),
    linewidth = 0.8
  ) +
  # Laplace PDF
  stat_function(
    aes(colour = "Laplace"),
    fun = function(x) (1/(2*bL))*exp(-abs(x-muL)/bL),
    linewidth = .8
  ) +
  scale_colour_manual(
    name = NULL,   
    values = c(
      "Normal"  = "#0072B2",  
      "Laplace" = "black"
    )
  ) +
  labs(
    x = expression(Delta*U[head]~"(m/s)"),
    y = "PDF"
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.title = element_text(size = 9 * 1.3),
    axis.text  = element_text(size = 9 * 1.3),
    legend.text = element_text(size = 9 * 1.3),
    legend.position = c(0.05, 0.95),  # 左上
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = "white",
      colour = NA
    )
  ) +
  xlim(-3, 3)



# Draw fig_3c QQplot
ggplot() +
  geom_point(
    data = qq_norm,
    aes(theo, samp, colour = "Normal", shape = "Normal"),
    size = 1
  ) +
  geom_point(
    data = qq_lap,
    aes(theo, samp, colour = "Laplace", shape = "Laplace"),
    size = 1
  ) +
  geom_abline(
    aes(colour = "Normal", linetype = "Normal"),
    intercept = 0, slope = 1,
    linewidth = 0.5
  ) +
  geom_abline(
    aes(colour = "Laplace", linetype = "Laplace"),
    intercept = 0, slope = 1,
    linewidth = 0.5
  ) +
  scale_colour_manual(
    name = NULL,  
    values = c("Normal" = "#0072B2", "Laplace" = "black")
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("Normal" = 1, "Laplace" = 4)
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("Normal" = "dashed", "Laplace" = "solid")
  ) +
  guides(
    colour   = guide_legend(override.aes = list(size = 2)),
    shape    = guide_legend(override.aes = list(size = 2)),
    linetype = guide_legend(override.aes = list(linewidth = 0.8))
  ) +
  coord_equal() +
  labs(
    x = "Theoretical quantile",
    y = "Empirical quantile"
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.title = element_text(size = 9 * 1.3),
    axis.text  = element_text(size = 9 * 1.3),
    legend.text = element_text(size = 9 * 1.3),
    legend.position = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = "white",
      colour = NA
    )
  )


# AIC Check

x_col <- "dU_1s"   # "
min_n <- 50          # 窓内の最小サンプル数（1Hzなら通常~600あるはず）

fit_norm_vs_laplace <- function(x){
  x <- x[is.finite(x)]
  n <- length(x)
  if(n < min_n) return(tibble(n=n, ll_norm=NA_real_, ll_lap=NA_real_,
                              deltaAIC=NA_real_, deltaLL_per=NA_real_))
  
  # --- Normal (MLE) ---
  muN <- mean(x)
  sdN <- sd(x)
  if(!is.finite(sdN) || sdN <= 0) return(tibble(n=n, ll_norm=NA_real_, ll_lap=NA_real_,
                                                deltaAIC=NA_real_, deltaLL_per=NA_real_))
  ll_norm <- sum(dnorm(x, mean = muN, sd = sdN, log = TRUE))
  
  # --- Laplace (MLE: median + mean abs dev about median) ---
  muL <- median(x)
  b   <- mean(abs(x - muL))
  if(!is.finite(b) || b <= 0) return(tibble(n=n, ll_norm=ll_norm, ll_lap=NA_real_,
                                            deltaAIC=NA_real_, deltaLL_per=NA_real_))
  ll_lap <- sum(-log(2*b) - abs(x - muL)/b)
  
  # パラメータ数は両方2（位置＋尺度）なので、AIC差は 2*(ll差) だけでOK
  # deltaAIC > 0 なら Laplace が良い
  deltaAIC    <- 2*(ll_lap - ll_norm)
  deltaLL_per <- (ll_lap - ll_norm)/n  # 1サンプルあたりの改善量
  
  tibble(n=n, ll_norm=ll_norm, ll_lap=ll_lap,
         deltaAIC=deltaAIC, deltaLL_per=deltaLL_per)
}

keep_win <- df_norm %>%
  group_by(block_10min) %>%
  summarise(mean_ws = mean(windspd, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(mean_ws >= 1) %>%
  pull(block_10min)

fit_df_all <- df_norm %>%
#  dplyr::filter(block_10min %in% keep_win) %>%
  group_by(block_10min) %>%
  group_modify(~ fit_norm_vs_laplace(.x[[x_col]])) %>%
  ungroup()

summary_fit <- fit_df_all %>%
  dplyr::filter(is.finite(deltaAIC)) %>%
  summarise(
    n_win = n(),
    prop_laplace_better = mean(deltaAIC > 0),
    deltaAIC_med = median(deltaAIC),
    deltaAIC_q25 = quantile(deltaAIC, 0.25),
    deltaAIC_q75 = quantile(deltaAIC, 0.75),
    deltaAIC_p95 = quantile(deltaAIC, 0.95),
    deltaLL_per_med = median(deltaLL_per)
  )

summary_fit

count_fit <- function(fit_df){
  fit_df %>%
    summarise(
      n_total = n(),
      n_finite = sum(is.finite(deltaAIC)),
      n_lap_win = sum(deltaAIC > 0, na.rm = TRUE),
      n_norm_win = sum(deltaAIC < 0, na.rm = TRUE),
      n_tie = sum(deltaAIC == 0, na.rm = TRUE),
      n_na = sum(!is.finite(deltaAIC))
    )
}


count_fit(fit_df_all)


ggplot(b_df, aes(x = b)) +
  geom_histogram(bins = 80) +
  labs(x = "b (m/s)", y = "Count") +
  theme_bw()
#################### Chapter 6 ##########################


# Caliculate fig_4
# ---- Set parameter ----
dt_list  <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)   # Δt [s]
win_sec  <- 600                              # window
fs       <- 10                               # 10 Hz sampling
min_mean_wind <- 1.0                         # minimum wind speed 
# 平均風速が 1m/s 未満の窓は，計測分解能およびノイズの影響が支配的となることを避けるため除外した．
min_pairs <- 1000                            # 



read_ulsa_csv <- function(file){
  dt <- fread(file)
  dt[, ts := as.POSIXct(ts, tz = "Asia/Tokyo")]
  setorder(dt, ts)
  dt <- dt[valid == 1]
  dt
}

add_10min_window <- function(dt, win_sec = 600){
  dt[, win_start := as.POSIXct(
    floor(as.numeric(ts) / win_sec) * win_sec,
    origin = "1970-01-01", tz = "Asia/Tokyo"
  )]
  dt
}

calc_increment <- function(x, lag_s, fs = 10){
  lag_n <- as.integer(round(lag_s * fs))   
  if (lag_n < 1 || length(x) <= lag_n) return(numeric(0))
  x[(lag_n + 1):length(x)] - x[1:(length(x) - lag_n)]
}

calc_b <- function(du){
  if (length(du) == 0) return(NA_real_)
  mu_hat <- median(du)
  mean(abs(du - mu_hat))
}

process_one_file <- function(file){
  
  dt <- read_ulsa_csv(file)
  dt <- add_10min_window(dt, win_sec)
  
  # ---- 窓ごとの平均風速を計算して、1m/s未満の窓を先に除外 ----
  win_mean <- dt[, .(mean_wind = mean(windspd_head, na.rm = TRUE)), by = win_start]
  keep_win <- win_mean[mean_wind >= min_mean_wind, win_start]
  
  dt <- dt[win_start %in% keep_win]
  if (nrow(dt) == 0) return(NULL)
  
  # ---- b(Δt) を計算（returnは使わず、空data.tableを返す） ----
  res <- dt[, {
    mean_wind <- mean(windspd_head, na.rm = TRUE)
    
    out <- rbindlist(lapply(dt_list, function(dtl){
      du <- calc_increment(windspd_head, dtl, fs)
      if (length(du) < min_pairs) return(NULL)
      
      data.table(
        delta_t   = dtl,
        N         = length(du),
        b         = calc_b(du),
        mean_wind = mean_wind
      )
    }), fill = TRUE)
    
    if (is.null(out) || nrow(out) == 0) {
      # data.tableのby計算では NULL を返さない方が安全
      out <- data.table(delta_t = numeric(0), N = integer(0), b = numeric(0), mean_wind = numeric(0))
    }
    out
  }, by = win_start]
  
  if (nrow(res) == 0) return(NULL)
  res[, file := basename(file)]
  res
}



# Load files in pole_ULSA_PRO/ULSA_PRO_10hz_20251207_20251218  Folder

folder <- file.path(dataset_root, "raw", "pole_ULSA_PRO",
                       "ULSA_PRO_10hz_20251207_20251218")
files <- list.files(folder, pattern="\\.csv$", full.names=TRUE)

result_all <- rbindlist(lapply(files, process_one_file), fill = TRUE)

plot_dt <- copy(result_all)
plot_dt[, delta_t_f := factor(delta_t, levels = sort(unique(delta_t)))]

# Show total windows  1window = 7Delta_t
nrow(plot_dt)/7
# 解析に際しては， 10 分窓ごとにデータを区切り，十分なサンプル数を確保しつつ， 
#∆t = 0.1–10 s の範囲における時間増分スケール解析が安定して行えるよう設定した（n = 347 窓）．


# Draw fig 4
#Box-and-whisker plots of the Laplace turbulence index b(∆t) for different time lags ∆t, evaluated over nonoverlapping 10-min windows. 
ggplot(plot_dt, aes(x = delta_t_f, y = b)) +
  geom_boxplot(
    width = 0.6
  ) +
  labs(
    x = expression(Delta*t*" (s)"),
    y = "b (m/s)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 9 * 1.3),
    axis.text  = element_text(size = 9 * 1.3),
    legend.text = element_text(size = 9 * 1.3)
  )






# （スケーリング指数）用設定 ----
dt_fit_range <- c(0.5, 10)   # フィットに使う Δt 範囲 [s]
min_points   <- 4            # 回帰に必要な最小Δt点数

estimate_zeta1 <- function(dt, b){
  
  ok <- is.finite(b) & b > 0
  dt <- dt[ok]
  b  <- b[ok]
  
  if (length(dt) < min_points) return(NA_real_)
  
  fit <- lm(log(b) ~ log(dt))
  coef(fit)[2]   # 傾き ζ1
}

estimate_zeta1_stats <- function(dt, b){
  
  ok <- is.finite(dt) & dt > 0 & is.finite(b) & b > 0
  dt <- dt[ok]
  b  <- b[ok]
  
  n_fit <- length(dt)
  if (n_fit < min_points) {
    return(list(
      zeta1 = NA_real_,
      r2    = NA_real_,
      r_fit = NA_real_,
      n_fit = n_fit
    ))
  }
  
  fit <- lm(log(b) ~ log(dt))
  s   <- summary(fit)
  
  zeta1 <- unname(coef(fit)[2])        # 傾き ζ1
  r2    <- unname(s$r.squared)         # 決定係数 R^2
  r_fit <- sign(zeta1) * sqrt(r2)      # log–log回帰の相関係数（符号付き）
  
  list(
    zeta1 = zeta1,
    r2    = r2,
    r_fit = r_fit,
    n_fit = n_fit
  )
}


zeta_dt <- result_all[
  delta_t >= dt_fit_range[1] &
    delta_t <= dt_fit_range[2] &
    N >= 1000,
  {
    st <- estimate_zeta1_stats(delta_t, b)
    .(
      zeta1     = st$zeta1,
      r2        = st$r2,       # （回帰の当てはまり）
      r_fit     = st$r_fit,    # （符号付き相関）
      n_fit     = st$n_fit,    # （回帰に使った点数）
      mean_wind = mean(mean_wind, na.rm = TRUE)
    )
  },
  by = .(file, win_start)
]

zeta_dt <- zeta_dt[is.finite(zeta1) & is.finite(r2)]

summary(zeta_dt$zeta1)
# 本データ（n = 347 窓）において推定された ζ1 は 0.3146–0.5000 の範囲に分布し，
#中央値は 0.4080，四分位範囲は 0.3848–0.4271 であった．

quantile(zeta_dt$zeta1, probs = c(0.25, 0.5, 0.75))



# b(1s) を抽出
b1_dt <- result_all[
  delta_t == 1 & N >= 1000,
  .(
    b1 = b
  ),
  by = .(file, win_start)
]

# ζ₁ と結合
zb_dt <- merge(
  zeta_dt,
  b1_dt,
  by = c("file", "win_start")
)

# 1) 図5全体：ζ1 と平均風速の相関（効果量）
pearson_res  <- cor.test(zeta_dt$mean_wind, zeta_dt$zeta1, method = "pearson")
spearman_res <- cor.test(zeta_dt$mean_wind, zeta_dt$zeta1, method = "spearman", exact = FALSE)

pearson_res
spearman_res

# ζ1 は平均風速と有意な負の相関を示した（Pearsonr = −0.474， 95%CI [−0.552, −0.388]， p = 7.53×10−21； Spearman ρ = −0.463， p = 7.97×10−20）．


cat("n =", nrow(zeta_dt), "\n")
cat("Pearson r =", unname(pearson_res$estimate),
    " 95%CI [", pearson_res$conf.int[1], ",", pearson_res$conf.int[2], "]",
    " p =", pearson_res$p.value, "\n")
cat("Spearman rho =", unname(spearman_res$estimate),
    " p =", spearman_res$p.value, "\n")

# 2) 各窓の log(b)~log(dt) 回帰の当てはまり（R^2）分布
summary(zeta_dt$r2)
quantile(zeta_dt$r2, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

# また，図 5 では各窓の回帰を表す指標として決定係数 R2 も併せて示した． R2 は 0.9166–0.9952（中央値 0.9733）と全体として非常に高い事が確認された．

cor.test(zeta_dt$mean_wind, zeta_dt$r2, method = "spearman", exact = FALSE)
cor.test(zeta_dt$mean_wind, zeta_dt$r2, method = "pearson")


#R^2と Wind speed の相関
zeta_dt[, wind_bin := cut(mean_wind, breaks = c(1,1.5,2,2.5,3,3.5,4), include.lowest=TRUE)]
zeta_dt[, .(
  n = .N,
  r2_med = median(r2, na.rm=TRUE),
  r2_q25 = quantile(r2, 0.25, na.rm=TRUE),
  r2_q75 = quantile(r2, 0.75, na.rm=TRUE)
), by = wind_bin]



# Draw fig5
# Relationship between the scaling exponent ζ1 and the mean wind speed evaluated over non-overlapping 10-min windows.

ggplot(zeta_dt, aes(mean_wind, zeta1, color=r2)) +
  geom_point(size=1.5) +
  labs(x="Wind speed (m/s)", y=expression(zeta[1]), color=expression(R^2)) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 9 * 1.3),
        axis.text  = element_text(size = 9 * 1.3),
        legend.text = element_text(size = 9 * 1.3))




### Chap 7


process_trisonica <- function(file,
                              fs = 10,
                              win_min = 10,
                              fc = 4.0,      # cutoff frequency [Hz]
                              order = 4) {
  
  dt <- fread(file)
  stopifnot(all(c("U", "V") %in% names(dt)))
  
  n_win <- fs * 60 * win_min
  n_all <- nrow(dt)
  
  if (n_all < n_win + fs) {
    warning(basename(file), " : データ長不足")
    return(NULL)
  }
  
  # ---- 中央10分 ----
  i_center <- floor(n_all / 2)
  i_start  <- i_center - floor(n_win / 2) + 1
  i_end    <- i_start + n_win - 1
  dt_win   <- dt[i_start:i_end]
  
  # ---- 処理前 WindSpeed ----
  dt_win[, WS_raw := sqrt(U^2 + V^2)]
  
  # ---- ローパス（U, V 個別）----
  dt_win[, U_f := lowpass_uv(U, fs, fc, order)]
  dt_win[, V_f := lowpass_uv(V, fs, fc, order)]
  
  # ---- 処理後 WindSpeed ----
  dt_win[, WS_filt := sqrt(U_f^2 + V_f^2)]
  
  # ---- 1秒中央値 ----
  dt_win[, time_sec := (seq_len(.N) - 1) %/% fs]
  
  uv_1s <- dt_win[
    ,
    .(
      U_1s = median(U_f, na.rm = TRUE),
      V_1s = median(V_f, na.rm = TRUE)
    ),
    by = time_sec
  ][order(time_sec)]
  
  # 共通：ラプラス尺度推定（中央値で中心化）
  laplace_b <- function(dx, min_n = 20){
    dx <- dx[is.finite(dx)]
    if(length(dx) < min_n) return(NA_real_)
    mu <- median(dx)
    mean(abs(dx - mu))
  }
  
  # ---- LTI（1秒代表値 → 1秒差分 → 中心化してb）----
  dx_u <- diff(uv_1s$U_1s)
  dx_v <- diff(uv_1s$V_1s)
  
  LTI_u <- laplace_b(dx_u)
  LTI_v <- laplace_b(dx_v)
  

  
  # ---- 時系列出力（比較・作図用）----
  ts_out <- dt_win[
    ,
    .(
      file,
      device = "Trisonica",
      time_sec = time_sec,
      WS_raw,
      WS_filt
    )
  ]
  
  list(
    summary = data.table(
      file  = basename(file),
      LTI_u = LTI_u,
      LTI_v = LTI_v,
      fc    = fc,
      order = order
    ),
    ts = ts_out
  )
}

lowpass_uv <- function(x, fs, fc, order = 4) {
  bf <- butter(order, fc / (fs / 2), type = "low")
  filtfilt(bf, x)
}



folder <- file.path(dataset_root, "raw", "drone_Trisonica")
files  <- list.files(folder, pattern="\\.csv$", full.names=TRUE)

res <- lapply(files, process_trisonica)

# ---- LTI（1ファイル=1値）----
Trisonica_LTI <- rbindlist(lapply(res, `[[`, "summary"), fill = TRUE)

# ---- 時系列（処理前・後）----
Trisonica_TS <- rbindlist(lapply(res, `[[`, "ts"), fill = TRUE)

#LTIを表示
Trisonica_LTI

#最も強い条件のデータを検索
max(Trisonica_LTI$LTI_u)
max(Trisonica_LTI$LTI_v)

#最も強い条件のデータを選択、表示
ggplot(Trisonica_TS[file == unique(file)[7]],
       aes(x = time_sec)) +
  geom_line(aes(y = WS_raw, color = "Raw"), alpha = 0.6) +
  geom_line(aes(y = WS_filt, color = "Filtered"), linewidth = 0.8) +
  labs(
    x = "Time [s]",
    y = "Wind speed [m/s]",
    color = "",
    title = "Effect of low-pass filtering on wind speed"
  ) +
  theme_bw()

# これらの 1 s 代表値を比較した結果， raw と filtered の間に顕著な差は認められず，
#ローパスフィルタ処理による影響が僅少であることが確認されたため，本節では両者の時系列を図として示すことは省略する．


#################

# ---- ファイルパス ----


file_drone <- file.path(dataset_root, "raw", "drone_Trisonica",
                        "wx_data_20251209_130549_130900-132859.csv")


file_pole <- file.path(dataset_root, "raw", "pole_ULSA_PRO",
                    "ULSA_PRO_10hz_20251207_20251218","wind_2025-12-09.csv")


# ---- 読み込み ----
df_drone <- read_csv(file_drone, show_col_types = FALSE)
df_pole  <- read_csv(file_pole,  show_col_types = FALSE)

df_drone <- df_drone %>%
  mutate(
    timestamp = force_tz(ymd_hms(timestamp), "Asia/Tokyo")
  )

df_pole <- df_pole %>%
  mutate(
    ts = force_tz(ymd_hms(ts), "Asia/Tokyo")
  )

df_drone <- as.data.frame(df_drone)
df_pole  <- as.data.frame(df_pole)




# ---- ドローンの中心時刻を基準に10分窓を作る ----
t_center <- df_drone$timestamp[ceiling(nrow(df_drone) / 2)]
t_start  <- t_center - minutes(5)
t_end    <- t_center + minutes(5)


# ---- 10分窓切り出し ----
df_drone_10min <- df_drone[
  df_drone$timestamp >= t_start &
    df_drone$timestamp <  t_end,
]

df_pole_10min <- df_pole[
  df_pole$ts >= t_start &
    df_pole$ts <  t_end,
]



fs <- 10      # Hz
fc <- 4       # Hz
bf <- butter(4, fc / (fs/2), type = "low")

# 10分窓の先頭を 1秒境界に合わせる（JST）
t0 <- floor_date(df_drone_10min$timestamp[1], "1 second")

df_drone_10min <- df_drone_10min %>%
  mutate(
    sec = (row_number() - 1) %/% fs,
    t1s = t0 + seconds(sec)
  )



df_drone_filt <- df_drone_10min %>%
  mutate(
    U = as.numeric(U),
    V = as.numeric(V),
    u_f = filtfilt(bf, U),
    v_f = filtfilt(bf, V)
  )

df_drone_10min <- df_drone_10min %>%
  mutate(
    U = as.numeric(U),
    V = as.numeric(V)
  )


# raw 1Hz：U,V の1秒中央値 → ws（sqrtは最後に1回）
drone_raw_1s <- df_drone_10min %>%
  group_by(t1s) %>%
  summarise(
    u_1s = median(U, na.rm = TRUE),
    v_1s = median(V, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ws   = sqrt(u_1s^2 + v_1s^2),
    type = "Drone raw"
  ) %>%
  select(t1s, ws, type)




drone_filt_1s <- df_drone_filt %>%
  group_by(t1s) %>%
  summarise(
    u_1s = median(u_f, na.rm = TRUE),
    v_1s = median(v_f, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ws   = sqrt(u_1s^2 + v_1s^2),
    type = "Drone filtered"
  ) %>%
  select(t1s, ws, type)

# pole 側も数値化
df_pole_10min <- df_pole_10min %>%
  mutate(
    windspd = as.numeric(windspd),
    sec = floor(as.numeric(difftime(ts, t0, units = "secs"))),
    t1s = t0 + seconds(sec)
  )

pole_1s <- df_pole_10min %>%
  group_by(t1s) %>%
  summarise(ws = median(windspd, na.rm = TRUE), .groups = "drop") %>%
  mutate(type = "Pole") %>%
  select(t1s, ws, type)


df_plot3 <- bind_rows(drone_raw_1s, drone_filt_1s, pole_1s)
df_plot2 <- bind_rows(drone_filt_1s, pole_1s)

# Fig6 7-3
ggplot(df_plot2, aes(x = t1s, y = ws, color = type)) +
  geom_line(linewidth = .5) +
  scale_x_datetime(timezone = "Asia/Tokyo") +
  scale_color_manual(
    values = c("Drone filtered"="#E69F00", "Pole"="#0072B2")
  ) +
  labs(x = "Time (JST)", y = "Wind speed (m/s)", color = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 9 * 1.3),
        axis.text  = element_text(size = 9 * 1.3),
        legend.text = element_text(size = 9 * 1.3))




############### 複数ファイルを処理する

dir_drone <- file.path(dataset_root, "raw", "drone_Trisonica")
dir_pole  <- file.path(dataset_root, "raw", "pole_ULSA_PRO","ULSA_PRO_10hz_select")

files_drone <- list.files(dir_drone, pattern="\\.csv$", full.names=TRUE)
files_pole  <- list.files(dir_pole,  pattern="\\.csv$", full.names=TRUE)

# 念のためソート（並び順を保証）
files_drone <- sort(files_drone)
files_pole  <- sort(files_pole)

process_one_pair <- function(file_drone, file_pole,
                             fs = 10, win_min = 10,
                             fc = 4, order = 4) {
  
  
  # --- 読み込み ---
  df_d <- read_csv(file_drone, show_col_types = FALSE) %>%
    mutate(timestamp = force_tz(ymd_hms(timestamp), "Asia/Tokyo")) %>%
    as.data.frame()
  
  df_p <- read_csv(file_pole, show_col_types = FALSE) %>%
    mutate(ts = force_tz(ymd_hms(ts), "Asia/Tokyo")) %>%
    as.data.frame()
  
  # --- 中央10分（行番号ベース）---
  n_win <- fs * 60 * win_min
  n_all <- nrow(df_d)
  if (n_all < n_win) return(NULL)
  
  i_center <- floor(n_all / 2)
  i_start  <- i_center - floor(n_win / 2) + 1
  i_end    <- i_start + n_win - 1
  
  df_d <- df_d[i_start:i_end, ]
  
  t_start <- df_d$timestamp[1]
  t_end   <- df_d$timestamp[nrow(df_d)]
  
  df_p <- df_p[df_p$ts >= t_start & df_p$ts <= t_end, ]
  
  # --- ローパス（Drone）---
  bf <- butter(order, fc / (fs/2), type="low")
  
  df_d <- df_d %>%
    mutate(
      U = as.numeric(U),
      V = as.numeric(V),
      u_f = filtfilt(bf, U),
      v_f = filtfilt(bf, V)
    )
  
  # --- 1 Hz 時刻（行番号ベース）---
  t0 <- floor_date(df_d$timestamp[1], "1 second")
  
  df_d <- df_d %>%
    mutate(
      sec = (row_number() - 1) %/% fs,
      t1s = t0 + sec
    )
  
  # --- Drone 1Hz ---
  drone_1s <- df_d %>%
    group_by(t1s) %>%
    summarise(
      u_1s = median(u_f, na.rm=TRUE),
      v_1s = median(v_f, na.rm=TRUE),
      .groups="drop"
    ) %>%
    mutate(
      ws = sqrt(u_1s^2 + v_1s^2),
      Device = "Drone"
    ) %>%
    select(ws, Device)
  
  # --- Pole 1Hz ---
  df_p <- df_p %>%
    mutate(
      windspd = as.numeric(windspd),
      sec = floor(as.numeric(difftime(ts, t0, units="secs"))),
      t1s = t0 + sec
    )
  
  # --- ドローン側の 1Hz 時刻（正） ---
  t1s_ref <- sort(unique(df_d$t1s))
  
  # --- Pole 側を「ドローンの t1s」に合わせる ---
  pole_1s <- lapply(t1s_ref, function(tt) {
    idx <- df_p$ts >= tt & df_p$ts < tt + 1
    if (!any(idx)) return(NULL)
    data.frame(
      t1s = tt,
      ws  = median(as.numeric(df_p$windspd[idx]), na.rm = TRUE)
    )
  }) |> 
    dplyr::bind_rows() |>
    mutate(Device = "Pole") |>
    select(ws, Device)
  
  bind_rows(drone_1s, pole_1s) %>%
    mutate(case = basename(file_drone))
}




df_all <- map2_dfr(
  files_drone,
  files_pole,
  process_one_pair
)



 # Pole基準で並べる
case_order <- df_all %>%
  dplyr::filter(Device == "Pole") %>%
  group_by(case) %>%
  summarise(mean_ws = mean(ws, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_ws) %>%
  pull(case)


case_labels <- tibble(case = case_order) %>%
  mutate(
    label = case %>%
      str_extract("\\d{8}_\\d{6}") %>%        # 20251209_130549
      ymd_hms(tz = "Asia/Tokyo") %>%
      format("%Y-%m-%d-%H-%M")
  )

df_plot_box <- df_all %>%
  mutate(
    case = factor(case, levels = case_order)
  )


# fig 7
ggplot(df_plot_box, aes(x = case, y = ws, fill = Device)) +
  geom_boxplot(
    position = position_dodge(width = 0.7),
    width = 0.6,
    linewidth = 0.4,
    outlier.size = 0.4
  ) +
  scale_x_discrete(
    labels = setNames(case_labels$label, case_labels$case)
  ) +
  scale_fill_manual(
    values = c("Drone" = "#E69F00", "Pole" = "#0072B2")
  ) +
  labs(
    x = "Start time",
    y = "Wind speed (m/s)",
    fill = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    ggplot(df_plot_box, aes(x = case, y = ws, fill = Device)) +
  geom_boxplot(
    position = position_dodge(width = 0.7),
    width = 0.6,
    linewidth = 0.4,
    outlier.size = 0.4
  ) +
  scale_x_discrete(
    labels = setNames(case_labels$label, case_labels$case)
  ) +
  scale_fill_manual(
    values = c("Drone" = "#E69F00", "Pole" = "#0072B2")
  ) +
  labs(
    x = "Start time",
    y = "Wind speed (m/s)",
    fill = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 9 * 1.3),
    axis.text  = element_text(size = 9 * 1.3),
    legend.text = element_text(size = 9 * 1.3)
  )
  )


############## bの比較　####################



process_one_pair <- function(file_drone, file_pole,
                             fs = 10, win_min = 10,
                             fc = 4, order = 4) {
  
  
  # =========================
  # b（Laplace Turbulence Index）計算関数
  # =========================
  calc_b_laplace <- function(x) {
    dx <- diff(x)
    dx <- dx[is.finite(dx)]
    if (length(dx) < 20) return(NA_real_)
    mean(abs(dx))
  }
  
  # =========================
  # 読み込み
  # =========================
  df_d <- read_csv(file_drone, show_col_types = FALSE) %>%
    mutate(timestamp = force_tz(ymd_hms(timestamp), "Asia/Tokyo")) %>%
    as.data.frame()
  
  df_p <- read_csv(file_pole, show_col_types = FALSE) %>%
    mutate(ts = force_tz(ymd_hms(ts), "Asia/Tokyo")) %>%
    as.data.frame()
  
  # =========================
  # 中央 10 分窓（Drone 基準）
  # =========================
  n_win <- fs * 60 * win_min
  n_all <- nrow(df_d)
  if (n_all < n_win) return(NULL)
  
  i_center <- floor(n_all / 2)
  i_start  <- i_center - floor(n_win / 2) + 1
  i_end    <- i_start + n_win - 1
  
  df_d <- df_d[i_start:i_end, ]
  
  t_start <- df_d$timestamp[1]
  t_end   <- df_d$timestamp[nrow(df_d)]
  
  df_p <- df_p[df_p$ts >= t_start & df_p$ts <= t_end, , drop = FALSE]
  
  
  # =========================
  # ローパスフィルタ（Drone）
  # =========================
  bf <- butter(order, fc / (fs / 2), type = "low")
  
  df_d <- df_d %>%
    mutate(
      U = as.numeric(U),
      V = as.numeric(V),
      u_f = filtfilt(bf, U),
      v_f = filtfilt(bf, V)
    )
  
  # =========================
  # Drone：1 Hz（中央値）
  # =========================
  t0 <- floor_date(df_d$timestamp[1], "1 second")
  
  df_d <- df_d %>%
    mutate(
      sec = (row_number() - 1) %/% fs,
      t1s = t0 + sec
    )
  
  drone_1s <- df_d %>%
    group_by(t1s) %>%
    summarise(
      u_1s = median(u_f, na.rm = TRUE),
      v_1s = median(v_f, na.rm = TRUE),
      .groups = "drop"
    )
  
  # =========================
  # Pole：Drone の t1s に合わせて 1 Hz（中央値）
  # =========================
  df_p <- df_p %>%
    mutate(
      windspd = as.numeric(windspd_head)
    )
  
  t1s_ref <- drone_1s$t1s
  
  pole_1s <- lapply(t1s_ref, function(tt) {
    idx <- df_p$ts >= tt & df_p$ts < tt + 1
    if (!any(idx)) return(NULL)
    data.frame(
      t1s = tt,
      ws  = median(df_p$windspd[idx], na.rm = TRUE)
    )
  }) %>%
    bind_rows()
  
  if (nrow(pole_1s) < 20) return(NULL)
  
  # =========================
  # b の計算（Δt = 1 s）
  # =========================
  b_du <- calc_b_laplace(drone_1s$u_1s)
  b_dv <- calc_b_laplace(drone_1s$v_1s)
  b_p  <- calc_b_laplace(pole_1s$ws)
  
  # =========================
  # 出力（U・V を分ける）
  # =========================
  tibble(
    case = basename(file_drone),
    b_pole = b_p,
    b_drone = c(b_du, b_dv),
    component = c("U", "V")
  )
}


df_b_compare <- map2_dfr(
  files_drone,
  files_pole,
  process_one_pair
)

case_map <- tibble(case = basename(files_drone)) %>%
  mutate(case_id = row_number())

df_b_compare <- df_b_compare %>%
  left_join(case_map, by = "case")

label_offset <- tibble::tibble(
  case_id = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
  dx = c( 0.03,  -0.03,  0.03,  0.03,  0.03,
          0.03,  0.03,  0.03,  0.03,  0.03, 0.03),
  dy = c(-0.02,  0, 0, 0, -.01,
         0.02, 0, -0.03, 0.025, -0.01, -0.03)
)

df_label <- df_b_compare %>%
  group_by(case_id) %>%
  summarise(
    base_x = mean(b_pole,  na.rm = TRUE),
    base_y = mean(b_drone, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(label_offset, by = "case_id") %>%
  mutate(
    label_x = base_x + dx,
    label_y = base_y + dy
  )




df_case_legend <- df_b_compare %>%
  distinct(case_id, case) %>%
  mutate(
    datetime = case %>%
      str_extract("\\d{8}_\\d{6}") %>%   # 例: 20250114_103000
      ymd_hms(tz = "Asia/Tokyo"),
    label = paste0(
      case_id, "  ",
      format(datetime, "%Y-%m-%d-%H")
    )
  ) %>%
  arrange(case_id)

legend_text <- paste(df_case_legend$label, collapse = "\n")

# Fig 8
ggplot(df_b_compare,
       aes(x = b_pole, y = b_drone, color = component)) +
  
  # 1:1 line（背面）
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  # U–V を結ぶ線（case 単位）
  geom_line(
    aes(group = case_id),
    color = "grey50",
    linewidth = 0.4,
    show.legend = FALSE
  ) +
  
  # 点
  geom_point(size = 3) +
  
  # ★ ペアごとに1つのラベル（黒）
  geom_label(
    data = df_label,
    aes(x = label_x, y = label_y, label = case_id),
    color = "black",
    size = 2.8,
    label.size = 0.1,
    fill = "white",
    inherit.aes = FALSE
  ) +
  
  # ★ 追加する番号対応凡例
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = legend_text,
    hjust = 1.05, vjust = -0.1,
    size = 3,
    color = "black"
  ) +
  
  coord_equal(clip = "off") +
  
  labs(
    x = expression(b[ pole ]~"(m/s)"),
    y = expression(b[ drone ]~"(m/s)"),
    color = "Component"
  ) +
  
  theme_bw(base_size = 13)








