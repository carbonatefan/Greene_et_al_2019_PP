%This script is published alongside Greene et al.,2019 'Early Cenozoic
%Decoupling of Climate and Carbonate Compensation Depth Trends'. Please
%cite usage accordingly.

%This script plots a subsidence curve for a single location and displays 
%the paleodepth for that location at a user-specified point in the past 
%using the subsidence equations from Cramer et al., 2009 'Ocean overturning
%since the Late Cretaceous: Inferences from a new benthic foraminiferal 
%isotope compliation', Paleoceanography 24(4), PA4216 with a simplified
%sediment unloading term.

clear all
close all

WaterDepth=4000; %Present day water depth in meters
BasementAge=100; %Basement age in millions of years
SedThickness=200; %Total sediment thickness in meters

sampleAge=0:1:BasementAge;
count=1;
 for t=1:length(sampleAge)
    if BasementAge >=20
        Id=WaterDepth-3051*(1-(8/pi^2)*exp(-0.0278*(BasementAge)))-0.66*(SedThickness);
		Pd_cr(count)=Id+3051*(1-(8/pi^2)*exp(-0.0278*(BasementAge-sampleAge(t))))+0.66*((SedThickness/BasementAge)*(BasementAge-sampleAge(t)));
    elseif BasementAge <20
        Id=WaterDepth-365*(BasementAge)^(1/2)-0.66*SedThickness;
        Pd_cr(count)=Id+365*sqrt(BasementAge-sampleAge(t))+0.66*((SedThickness/BasementAge)*(BasementAge-sampleAge(t)));
    end
    count=count+1;
 end

PaleoDepth=find(sampleAge == 50); %Enter age (in millions of years)for paleodepth reconstruction
disp(Pd_cr(PaleoDepth))

plot(sampleAge,Pd_cr)
set(gca,'Xdir','reverse');
set(gca,'Ydir','reverse');
xlabel('Age (Ma)');
ylabel('Paleodepth (m)');

