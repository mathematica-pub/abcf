if (!exists('overall')) {
  overall <- readRDS('Data/all-overall.RDS')
}
if (!exists('mixing')) {
  mixing <- readRDS('Data/all-mixing.RDS')
}
if (!exists('ex_smy')) {
  ex_smy <- readRDS('Data/all-exemplar-summy.RDS')
}
if (!exists('msc_list')) {
  msc_list <- readRDS('Data/all-misc.RDS')
}
if (!exists('calib')) {
  calib <- readRDS('Data/all-calib.RDS')
}
if (!exists('scenarios')) {
  scenarios <- readRDS('Data/_scenarios.RDS')
}