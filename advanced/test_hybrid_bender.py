from shadow4.beamline.s4_beamline import S4Beamline
from shadow4.tools.logger import is_verbose, is_debug, set_verbose

#set_verbose()

beamline = S4Beamline()

# electron beam
from shadow4.sources.s4_electron_beam import S4ElectronBeam
electron_beam = S4ElectronBeam(energy_in_GeV=6,energy_spread=0.00138,current=0.2)
electron_beam.set_sigmas_all(sigma_x=1.48e-05,sigma_y=2.8e-06,sigma_xp=3.7e-06,sigma_yp=1.5e-06)

# magnetic structure
from shadow4.sources.undulator.s4_undulator_gaussian import S4UndulatorGaussian
source = S4UndulatorGaussian(
    period_length     = 0.02,     # syned Undulator parameter (length in m)
    number_of_periods = 235.2, # syned Undulator parameter
    photon_energy     = 8800.0, # Photon energy (in eV)
    delta_e           = 3.0, # Photon energy width (in eV)
    ng_e              = 100, # Photon energy scan number of points
    flag_emittance    = 1, # when sampling rays: Use emittance (0=No, 1=Yes)
    flag_energy_spread = 0, # when sampling rays: Use e- energy spread (0=No, 1=Yes)
    harmonic_number    = 1, # harmonic number
    flag_autoset_flux_central_cone  = 0, # value to set the flux peak
    flux_central_cone  = 10000000000.0, # value to set the flux peak
    )


# light source
from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource
light_source = S4UndulatorGaussianLightSource(name='GaussianUndulator', electron_beam=electron_beam, magnetic_structure=source,nrays=50000,seed=5676561)
beam = light_source.get_beam()

beamline.set_light_source(light_source)

# optical element number XX
from syned.beamline.shape import Rectangle
boundary_shape = Rectangle(x_left=-0.001, x_right=0.001, y_bottom=-0.0005, y_top=0.0005)

from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
optical_element = S4Screen(name='FE Mask', boundary_shape=boundary_shape,
    i_abs=0, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
    i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

from syned.beamline.element_coordinates import ElementCoordinates
coordinates = ElementCoordinates(p=25.5, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.14159)
from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

beam, footprint = beamline_element.trace_beam()

beamline.append_beamline_element(beamline_element)

# optical element number XX
from syned.beamline.shape import Rectangle
boundary_shape = Rectangle(x_left=-0.0005, x_right=0.0005, y_bottom=-0.0005, y_top=0.0005)

from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
optical_element = S4Screen(name='WB Slits', boundary_shape=boundary_shape,
    i_abs=0, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
    i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

from syned.beamline.element_coordinates import ElementCoordinates
coordinates = ElementCoordinates(p=1.7, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.14159)
from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

beam, footprint = beamline_element.trace_beam()

beamline.append_beamline_element(beamline_element)

# optical element number XX
from syned.beamline.shape import Rectangle
boundary_shape = Rectangle(x_left=-0.005, x_right=0.005, y_bottom=-0.17, y_top=0.17)
   
from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror
optical_element = S4PlaneMirror(name='Plane Mirror',boundary_shape=boundary_shape,
    f_reflec=0,f_refl=0,file_refl='<none>',refraction_index=0.99999+0.001j,
    coating_material='Si',coating_density=2.33,coating_roughness=0)

from syned.beamline.element_coordinates import ElementCoordinates
coordinates = ElementCoordinates(p=0.8, q=0, angle_radial=1.5678, angle_azimuthal=4.71239, angle_radial_out=1.5678)
movements = None
from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirrorElement
beamline_element = S4PlaneMirrorElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)

beam, mirr = beamline_element.trace_beam()

beamline.append_beamline_element(beamline_element)

# optical element number XX
from syned.beamline.shape import Rectangle
boundary_shape = Rectangle(x_left=-0.005, x_right=0.005, y_bottom=-0.17, y_top=0.17)
   
from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror
optical_element = S4PlaneMirror(name='Plane Mirror',boundary_shape=boundary_shape,
    f_reflec=0,f_refl=0,file_refl='<none>',refraction_index=0.99999+0.001j,
    coating_material='Si',coating_density=2.33,coating_roughness=0)

from syned.beamline.element_coordinates import ElementCoordinates
coordinates = ElementCoordinates(p=2.6, q=0, angle_radial=1.5678, angle_azimuthal=3.14159, angle_radial_out=1.5678)
movements = None
from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirrorElement
beamline_element = S4PlaneMirrorElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)

beam, mirr = beamline_element.trace_beam()

beamline.append_beamline_element(beamline_element)

# optical element number XX
from syned.beamline.shape import Rectangle
boundary_shape = Rectangle(x_left=-0.0005, x_right=0.0005, y_bottom=-0.0005, y_top=0.0005)

from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
optical_element = S4Screen(name='PB Slits', boundary_shape=boundary_shape,
    i_abs=0, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
    i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

from syned.beamline.element_coordinates import ElementCoordinates
coordinates = ElementCoordinates(p=2.3, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.14159)
from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

beam, footprint = beamline_element.trace_beam()

beamline.append_beamline_element(beamline_element)

# optical element number XX
from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
optical_element = S4PlaneCrystal(name='Plane Crystal',
    boundary_shape=None, material='Si',
    miller_index_h=1, miller_index_k=1, miller_index_l=1,
    f_bragg_a=False, asymmetry_angle=0.0,
    is_thick=1, thickness=0.001,
    f_central=1, f_phot_cent=0, phot_cent=8800.0,
    file_refl='bragg.dat',
    f_ext=0,
    material_constants_library_flag=0, # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
    )
from syned.beamline.element_coordinates import ElementCoordinates
coordinates = ElementCoordinates(p=0.8, q=0.01, angle_radial=1.34416, angle_azimuthal=3.14159, angle_radial_out=1.34416)
movements = None
from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
beamline_element = S4PlaneCrystalElement(optical_element=optical_element,coordinates=coordinates, movements=movements, input_beam=beam)

beam, mirr = beamline_element.trace_beam()

beamline.append_beamline_element(beamline_element)

# optical element number XX
from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
optical_element = S4PlaneCrystal(name='Plane Crystal',
    boundary_shape=None, material='Si',
    miller_index_h=1, miller_index_k=1, miller_index_l=1,
    f_bragg_a=False, asymmetry_angle=0.0,
    is_thick=1, thickness=0.001,
    f_central=1, f_phot_cent=0, phot_cent=8800.0,
    file_refl='bragg.dat',
    f_ext=0,
    material_constants_library_flag=0, # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
    )
from syned.beamline.element_coordinates import ElementCoordinates
coordinates = ElementCoordinates(p=0, q=0.99, angle_radial=1.34416, angle_azimuthal=3.14159, angle_radial_out=1.34416)
movements = None
from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
beamline_element = S4PlaneCrystalElement(optical_element=optical_element,coordinates=coordinates, movements=movements, input_beam=beam)

beam, mirr = beamline_element.trace_beam()

beamline.append_beamline_element(beamline_element)

# optical element number XX

from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
optical_element = S4Empty(name='to Mono Slits')

from syned.beamline.element_coordinates import ElementCoordinates
coordinates = ElementCoordinates(p=20, q=1, angle_radial=0, angle_azimuthal=4.71239, angle_radial_out=3.14159)
from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

beam, mirr = beamline_element.trace_beam()

beamline.append_beamline_element(beamline_element)

# optical element number XX
from syned.beamline.shape import Rectangle
boundary_shape = Rectangle(x_left=-0.0121, x_right=0.0121, y_bottom=-0.15, y_top=0.15)
        
from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
optical_element = S4EllipsoidMirror(name='Flexural Hinge Bender Mirror', boundary_shape=boundary_shape,
    surface_calculation=0,
    min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
    p_focus=55.700000, q_focus=0.640000, grazing_angle=0.002500,
    is_cylinder=1, cylinder_direction=0, convexity=1,
    f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999+0.001j,
    coating_material='Si', coating_density=2.33, coating_roughness=0)


from shadow4_advanced.benders.s4_flexural_hinge_bender_ellipsoid_mirror import S4FlexuralHingeBenderEllipsoidMirror

ellipsoid_mirror = optical_element
from srxraylib.profiles.benders.bender_io import BenderMovement

fit_to_focus_parameters = None
bender_movement = BenderMovement(position_upstream=168.0, position_downstream=252.0)
M1              = 0.6659
e               = 0.81404
ratio           = 0.8911247935125394
from srxraylib.profiles.benders.flexural_hinge_bender_manager import CalibrationParameters

calibration_parameters = CalibrationParameters(parameters_upstream=[0.009, 0.05], parameters_downstream=[0.006, 0.05])

optical_element = S4FlexuralHingeBenderEllipsoidMirror(ellipsoid_mirror=ellipsoid_mirror,
                                                       figure_error_data_file='KBV.hdf5',
                                                       bender_bin_x=10,
                                                       bender_bin_y=500,
                                                       E=131000000000.0,
                                                       h=0.0095,
                                                       M1=M1,
                                                       e=e,
                                                       ratio=ratio,
                                                       bender_shape=0,
                                                       bender_type=1,
                                                       calibration_parameters=calibration_parameters,
                                                       fit_to_focus_parameters=fit_to_focus_parameters,
                                                       bender_movement=bender_movement)

bender_data = optical_element.get_bender_data()
from syned.beamline.element_coordinates import ElementCoordinates
coordinates = ElementCoordinates(p=0, q=0.26, angle_radial=1.5683, angle_azimuthal=0, angle_radial_out=1.5683)
movements = None
from shadow4_advanced.benders.s4_flexural_hinge_bender_ellipsoid_mirror import S4FlexuralHingeBenderEllipsoidMirrorElement
beamline_element = S4FlexuralHingeBenderEllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)

beam, mirr = beamline_element.trace_beam()

beamline.append_beamline_element(beamline_element)

# optical element number XX
from syned.beamline.shape import Rectangle
boundary_shape = Rectangle(x_left=-0.01185, x_right=0.01185, y_bottom=-0.08, y_top=0.08)
        
from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
optical_element = S4EllipsoidMirror(name='Flexural Hinge Bender Mirror', boundary_shape=boundary_shape,
    surface_calculation=0,
    min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
    p_focus=55.960000, q_focus=0.380000, grazing_angle=0.002500,
    is_cylinder=1, cylinder_direction=0, convexity=1,
    f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999+0.001j,
    coating_material='Si', coating_density=2.33, coating_roughness=0)


from shadow4_advanced.benders.s4_flexural_hinge_bender_ellipsoid_mirror import S4FlexuralHingeBenderEllipsoidMirror

ellipsoid_mirror = optical_element
from srxraylib.profiles.benders.bender_io import BenderMovement

fit_to_focus_parameters = None
bender_movement = BenderMovement(position_upstream=103.5, position_downstream=135.8)
M1              = 0.7566
e               = 0.83122
ratio           = 0.9653713983610891
from srxraylib.profiles.benders.flexural_hinge_bender_manager import CalibrationParameters

calibration_parameters = CalibrationParameters(parameters_upstream=[0.025, 0.05], parameters_downstream=[0.019, 0.05])

optical_element = S4FlexuralHingeBenderEllipsoidMirror(ellipsoid_mirror=ellipsoid_mirror,
                                                       figure_error_data_file='KBH.hdf5',
                                                       bender_bin_x=10,
                                                       bender_bin_y=500,
                                                       E=131000000000.0,
                                                       h=0.0085,
                                                       M1=M1,
                                                       e=e,
                                                       ratio=ratio,
                                                       bender_shape=0,
                                                       bender_type=1,
                                                       calibration_parameters=calibration_parameters,
                                                       fit_to_focus_parameters=fit_to_focus_parameters,
                                                       bender_movement=bender_movement)

bender_data = optical_element.get_bender_data()
from syned.beamline.element_coordinates import ElementCoordinates
coordinates = ElementCoordinates(p=0, q=0.38, angle_radial=1.5683, angle_azimuthal=4.71239, angle_radial_out=1.5683)
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=1e-05, rotation_x=0, rotation_y=0, rotation_z=0)
from shadow4_advanced.benders.s4_flexural_hinge_bender_ellipsoid_mirror import S4FlexuralHingeBenderEllipsoidMirrorElement
beamline_element = S4FlexuralHingeBenderEllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)

beam, mirr = beamline_element.trace_beam()

beamline.append_beamline_element(beamline_element)

# optical element number XX
from hybrid_methods.coherence.hybrid_screen import StdIOHybridListener
from shadow4_advanced.hybrid.s4_hybrid_screen import S4HybridBeam, S4HybridOE, HybridInputParameters, S4HybridScreenElement
from shadow4_advanced.hybrid.s4_hybrid_screen import S4HybridScreen

calculation_type=7
hybrid_screen = S4HybridScreen(calculation_type)


additional_parameters = {}
kb_mirror_1_element = beamline.get_beamline_element_at(-2)
kb_mirror_2_element = beamline.get_beamline_element_at(-1)
input_parameters = HybridInputParameters(listener=StdIOHybridListener(),
                                         beam=S4HybridBeam(beam=[kb_mirror_2_element.get_input_beam(), beam]),
                                         optical_element=S4HybridOE(optical_element=[kb_mirror_1_element, kb_mirror_2_element]),
                                         diffraction_plane=3,
                                         propagation_type=1,
                                         n_bins_x=100,
                                         n_bins_z=100,
                                         n_peaks=50,
                                         fft_n_pts=100000,
                                         analyze_geometry=False,
                                         random_seed=None,
                                         **additional_parameters)

beamline_element = S4HybridScreenElement(hybrid_screen=hybrid_screen, hybrid_input_parameters=input_parameters)

beam, mirr, _ = beamline_element.trace_beam()

beamline.append_beamline_element(beamline_element)


# test plot
if True:
   from srxraylib.plot.gol import plot_scatter
   plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1), title='(Intensity,Photon Energy)', plot_histograms=0)
   plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1), title='(X,Z) in microns')
