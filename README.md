# FreThaw1D

Freezing and thawing of soils affect a wide range of biogeochemical and hydrological ([Walvoord and Kurylyk, 2016)](https://doi.org/10.2136/vzj2016.01.0010); [Schuur
et al., 2015](https://doi.org/10.1038/nature14338)) processes and interact with engineered structures in cold regions. This foundational importance together with the vast extent of cold regions globally makes the simulation of freezing and thawing soil an important and well-researched topic ([Streletskiy et al., 2019](https://doi.org/10.1088/1748-9326/aaf5e6); [Walvoord and Kurylyk, 2016)](https://doi.org/10.2136/vzj2016.01.0010); [Harris et al., 2009)](https://doi.org/10.1016/j.earscirev.2008.12.002)). The interest in understanding and anticipating the consequences of permafrost thaw driven by anthropogenic climate change has added additional urgency.

Studying frozen soil can be done at different timescales. At short timescale, the latent heat of fusion of water determines a time-lag in the cooling and warming of the soil thus affecting the surface heat fluxes. Studies have shown that proper frozen soil scheme help improve land surface and climate model simulation ([Viterbo et al., 1999](https://doi.org/10.1002/qj.49712555904); [Smirnova et al., 2000](https://doi.org/10.1029/1999JD901047)). 
On the other hand, permafrost evolution is characterized by a long timescale and its evolution is manly controlled by climate change. Permafrost degradation has significant impact on the ecosystems and on the social and economic life in cold regions ([Bao et al., 2016](https://doi.org/10.1002/2015JD024451))
Therefore, it is desirable to have a numerical model that can efficiently simulate frozen soil at both short timescale and long timescale.
At present, a common problem of physical-based frozen soil model concerns the treatment of the nonlinear behavior of the energy and enthalpy as function of temperature, water content and other governing variables. Because of this nonlinearity, the model may fail to compute the correct solution, or a short time step must be chosen leading to an increase of the computational cost.

`FreThaw1D` model solves the so called enthalpy - or conservative - formulation. The novelty regards the use of the nested Newton-Casulli-Zanolli (NCZ) algorithm, ([Casulli and Zanolli](https://doi.org/10.1137/100786320)) to linearize the nonlinear system resulting from the approximation of the governing equation. Using the enthalpy formulation and the NCZ algorithm the convergence of the solver is guarateed for any time step size. This is a key feature since the integration time step can be choosen accordingly to the time scale of the processes to study, from seconds to days.


![Alt text](docs/frozensoil.jpg?raw=true "Title")

#Acknowledgements

This work has been partially supported by a Ph.D. grant by the Department of Civil, Environmental and Mechanics Engineering at the University of Trento and by the NSERC PermafrostNet project. 
