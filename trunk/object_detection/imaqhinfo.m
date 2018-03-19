imaqhwinfo
dev_info=imaqhwinfo('winvideo',1)
%info=imaqhwinfo('winvideo')
   celldisp(dev_info.SupportedFormats)
  vid=videoinput('winvideo',1,'YUY2_320x240');
%Open Figure
hFigure=figure(1);
%set parameters for video
%Acquire only one frame each time
triggerconfig(vid,'Manual');
set(vid,'framespertrigger',1);
%Go on forever untill stopped
set(vid,'triggerrepeat',Inf);
%Get a grayscale image
set(vid,'ReturnedColorSpace','grayscale');
start(vid);
preview(vid);
for i=1:5
    trigger(vid);
    im=getdata(vid);
    figure,imshow(im(:,:,:,i));
end
stop(vid),delete(vid),clear(vid)