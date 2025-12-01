from equipment import Telescope, Camera, Optics

# Telescope: 80mm aperture, 600mm focal length, with 0.63x reducer
t = Telescope("ED80", 80, 600)
c = Camera("DSLR", 22, 14, 4000, 2600)
opt = Optics(t, c, focal_reducer=0.63)

print("Effective Focal Length (mm):", opt.effective_focal_length_mm)
print("FOV (deg):", opt.fov_deg)
print("Pixel scale (arcsec/pixel):", opt.pixel_scale_arcsec)
