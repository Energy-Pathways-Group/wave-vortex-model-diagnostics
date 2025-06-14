## Wave-Vortex Diagnostics

The `WVDiagnostics` is a class for quickly analyzing output from the wave-vortex model, with a focus on examining the energy fluxes from forcing and triads.

### Create a diagnostics file

The diagnostics class can create some figures without a diagnostics file, but in general you will want to create a diagnostics file first.

You instantiate a `WVDiagnostics` by pointing to model output,
```
wvd = WVDiagnostics(pathToWVModelOutput);
```
and then create a diagnostics file,
```
wvd.createDiagnostics();
```

There are several options for creating diagnostics file, including setting a stride (skipping time points), and explicitly computing the effect of antialising.

### Core energy flux functions

There are really two core energy flux functions upon which all the other energy flux functions are built:
```
forcing_fluxes = forcingFluxes(self,options);
inertial_fluxes = inertialFluxes(self,options);
```
both of which return arrays of structs. The forcing fluxes will have one entry for each forcing in the model, and the inertial fluxes will have one entry for each triad (explained below). The structs always have a `name` field, `fancyName` name field, and a field for each energy reservoir (explained below). The fields for each energy reservoir contain all the data, usually a matrix of dimension `[j kRadial t]`.

*Energy reservoirs* are defined by the `EnergyReservoir` class, and encompass several common combinationgs of energy reservoirs which might be useful for analysis. The two core fluxes functions have default energy reservoirs, but you can override this, e.g.,
```
energyReservoirs = [EnergyReservoir.geostrophic_kinetic, EnergyReservoir.geostrophic_potential, EnergyReservoir.wave, EnergyReservoir.total];
forcing_fluxes = wvd.forcingFluxes(energyReservoirs = energyReservoirs);
inertial_fluxes = wvd.inertialFluxes(energyReservoirs = energyReservoirs);
```
will separate out the geostrophic energy into kinetic energy and potential energy.

*Triad components* are very similar to energy reservoirs. By default the model has four primary components: mda, inertial, internal gravity wave, and geostrophic. However, this creates a mess of triads so the `TriadFlowComponent` class (part of this diagnostics tool) lets you optionally group together some of the flow components. For example,
```
forcing_fluxes = wvd.forcingFluxes(triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]);
```
combines the geostrophic and mda components together, and also combines the igw and io components together.
