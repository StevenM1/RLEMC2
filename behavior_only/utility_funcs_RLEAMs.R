## functions
s_left <- function(d) {
  gsub("([[:alpha:]]+)-([0-9]+\\.[0-9]+)_([[:alpha:]]+)-([0-9]+\\.[0-9]+)", '\\1', as.character(d$S))
}
p_left <- function(d) {
  as.numeric(gsub("([[:alpha:]]+)-([0-9]+\\.[0-9]+)_([[:alpha:]]+)-([0-9]+\\.[0-9]+)", '\\2', as.character(d$S)))
}
s_right <- function(d) {
  gsub("([[:alpha:]]+)-([0-9]+\\.[0-9]+)_([[:alpha:]]+)-([0-9]+\\.[0-9]+)", '\\3', as.character(d$S))
}
p_right <- function(d) {
  as.numeric(gsub("([[:alpha:]]+)-([0-9]+\\.[0-9]+)_([[:alpha:]]+)-([0-9]+\\.[0-9]+)", '\\4', as.character(d$S)))
}
p_low <- function(d) {
  pmin(d$p_left, d$p_right)
}
p_high <- function(d) {
  pmax(d$p_left, d$p_right)
}
s_low <- function(d) {
  ifelse(d$p_low==d$p_left, d$s_left, d$s_right)
}
s_high <- function(d) {
  ifelse(d$p_high==d$p_left, d$s_left, d$s_right)
}
lS <- function(d) {
  ## latent STIMULUS that matches an accumulator
  factor(d[cbind(1:nrow(d), match(paste0('s_', d$lR), colnames(d)))])
}
lRS <- function(d) {
  ## latent STIMULUS that was chosen
  ifelse(d$R=='left', d$s_left, d$s_right)
}
## latent DIRECTION that matches an accumulator
correct_direction <- function(d) {
  factor(ifelse(d$p_right>d$p_left, 'right', 'left'), levels=levels(d$R))
}
