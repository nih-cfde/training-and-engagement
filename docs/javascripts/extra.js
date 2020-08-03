// Turn any image with class attribute of "zoom-default" to a medium-zoom image
img = document.getElementsByTagName('img');
for(i = 0; i < img.length; i++) {
    img[i].className += " zoom-default";
}
const zoomDefault = mediumZoom('.zoom-default', {margin: 20,
  scrollOffset: 0,
  background: '#00000066'})
