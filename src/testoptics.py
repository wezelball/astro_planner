from optics import Optics

opt = Optics(
    telescope_aperture_mm=150,
    telescope_focal_length_mm=1000,
    camera_sensor_width_mm=22.3,
    camera_sensor_height_mm=14.9,
    camera_pixel_width=6000,
    camera_pixel_height=4000,
    focal_reducer=0.63
)

print("FOV (deg):", opt.fov_deg())
print("Pixel scale (arcsec/pixel):", opt.pixel_scale_arcsec())
print("FOV fill for 1.5 deg object:", opt.fov_fill_fraction(1.5))
