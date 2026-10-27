# RubinsForge

Rubin's repeated views of the sky are a gift to time-domain astronomy and a challenge for anything that leaves a trail across an exposure. RubinsForge brings together our LSST image simulations and satellite-trail analyses, from modeled motion to contamination in survey images.

Our `paths.py` module reads `RUBINSFORGE_PATH`, which defaults to this checkout. It exposes local `data/` and `visuals/` directories. Git ignores both.

The path module uses TDpy's shared repository-path utility. Install its pinned version before running the notebooks or scripts.

```bash
python -m pip install "tdpy @ git+https://github.com/tansudaylan/tdpy.git@35df1703c6b63d6a62f7d44b12cf5928b5d355ed"
```

For OpSim notebooks, `RUBIN_SIM_DATA_DIR` can point to a Rubin simulation data tree elsewhere. If it is unset, the notebooks look in `data/rubin_sim_data/`.
