classdef EnergyReservoir
   properties
      hasMDA logical
      hasGeostrophicKinetic logical
      hasGeostrophicPotential logical
      hasInertial logical
      hasInternalGravityWave logical
      name
      fancyName
   end
   properties (Dependent)
        vectorContents
   end
   methods
       function er = EnergyReservoir(hasMDA, hasGeostrophicKinetic, hasGeostrophicPotential, hasInertial, hasInternalGravityWave,name,fancyName)
         er.hasMDA = hasMDA;
         er.hasGeostrophicKinetic = hasGeostrophicKinetic;
         er.hasGeostrophicPotential = hasGeostrophicPotential;
         er.hasInertial = hasInertial;
         er.hasInternalGravityWave = hasInternalGravityWave;
         er.name = name;
         er.fancyName = fancyName;
       end

       function v = get.vectorContents(self)
            v = [self.hasMDA,self.hasGeostrophicKinetic,self.hasGeostrophicPotential,self.hasInertial,self.hasInternalGravityWave].';
       end
   end
   enumeration
      mda                       (true,false,false,false,false,"te_mda","mean density anomaly")
      geostrophic_kinetic       (false,true,false,false,false,"ke_g","geostrophic kinetic")
      geostrophic_potential     (false,false,true,false,false,"pe_g","geostrophic potential")
      geostrophic_potential_mda (true,false,true,false,false,"pe_g","geostrophic potential + mda")
      geostrophic               (false,true,true,false,false,"te_g","geostrophic")
      geostrophic_mda           (true,true,true,false,false,"te_gmda","geostrophic + mda")
      igw                       (false,false,false,false,true,"te_igw","internal gravity wave")
      io                        (false,false,false,true,false,"te_io","inertial")
      wave                      (false,false,false,true,true,"te_wave","wave")
      total                     (true,true,true,true,true,"te_quadratic","total quadratic")
   end
end