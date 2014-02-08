starcalib
=========

Lens calibration program using the night sky as calibration target. Useful for calibrating
wide-angle lenses, particularly when it is undesireable to refocus the lens after
calibration, which precludes the use of indoor calibration targets.

Usage: starcalib [OPTION] file.pnm

Calibrates a camera and optics by matching an image of the night sky against a star
catalogue

Options:
  -h, --help         Print this help
  -f, --foclen=MM    Focal length in mm (required)
  -p, --pitch=UM     Detector pixel pitch in microns (-p or -w required)
  -w, --width=MM     Effective detector width in mm (-p or -w required)
  -n, --nparam=n     Number of parameters to use in model, must be >= 3
                     (default = 11)
  -t, --threshold=T  Detection threshold, as ratio of star pixels to all pixels
                     (default = 5e-4)
  -o, --output=F     Name of output image file
Output:
  xo yo r rx ry r3 r3x r3y r5 r5x r5y ...

Corrected pixel positions (xc, yc) can be calculated from measued positions (xm, ym)
using the calibration parameters. Coordinates should be translated such that (0,0)
is at the center of the image and scaled so that the increment from one pixel to the
next is pitch/foclen (see options). Positive y-axis points upwards. Sample code:

  X = xm + xo;
  Y = ym + yo;
  R2 = X*X + Y*Y;
  c = r + X * rx + Y * ry + R2 * r3 + R2*X * r3x + R2*Y * r3y + R2*R2 * r5 ...
  xc = xm + X*c;
  yc = ym + Y*c;
  
