# R script for image processing from frames extracted from a video file from a chemical reaction
# inside the directory "/media/run"
library("imager")

## Define ROI dimensions horizontal (x) and vertical (y)
plot(load.image(file="/media/run/00000001.jpg"))
plot(imsub(load.image("/media/run/00000001.jpg"), x<560, x>160, y<800, y>400))
# load all image files
file_list <- list.files(path="/media/run/", full.names = TRUE)
# extract the mean value of the S channel
frames.mean  = list()
for (i in file_list) {
  frame.x = imsub(load.image(i), x<560, x>160, y<800, y>400)
  s.mean = mean((RGBtoHSV(frame.x))[,,2])
  frames.mean[[i]] = c(s.mean)
}
run <- frames.mean
# save the S mean
save(run,file = "Smean.RData")     
