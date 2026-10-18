"""Post-processing: analysis, plots, movies and interactive viewer."""
from .analysis import Run, rdf, msd, velocity_distribution, density_profile, radial_profile, coupling, fall_times  # noqa: F401


def __getattr__(name):  # lazy: render needs PyVista, which is optional
    if name == "render":
        import importlib
        return importlib.import_module(__name__ + ".render")
    raise AttributeError(name)
