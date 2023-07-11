## BayesTraitsV4-compatible parse_schedule from btwV2
## Kate Snyder
## 7/5/2023


parse_scheduleKTS <- function (file) {
  out <- scan(file = file, what = "c", quiet = T, sep = "\n")
  header_start <- grep("Rate Tried", out)
  #schedule <- read.table(file, nrow = (header_start - 1)) # , sep = "\t", header = FALSE
  #names(schedule) <- c("operator", "percent_tried")
  schedule = readLines(file, n = (header_start - 1))
  header <- read.table(file, skip = (header_start - 1), sep = "\t", 
                       header = TRUE)
  header <- header[, -ncol(header)]
  Schedule <- list(schedule = schedule, header = header)
  return(Schedule)
}
