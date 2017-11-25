/*
 * Copyright (C) 2017 user1
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package psu.ets.psubams.model;

/**
 *
 * @author user1
 */
/*
 * Copyright (C) 2015 Toby N. Carlson
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

import uk.ac.aber.simsphere.*;
import java.io.*;
import com.csvreader.CsvWriter;
import static psu.ets.psubams.model.Constants.CONVOLUTION_COUNT;
import static psu.ets.psubams.model.Constants.OUTPUT_COUNT;
import static psu.ets.psubams.model.Constants.SOUNDING_COUNT;
import static psu.ets.psubams.model.Constants.SOUNDING_LEVELS;
import static psu.ets.psubams.model.SeriesNames.OUT_AIR_TEMP_10M;
import static psu.ets.psubams.model.SeriesNames.OUT_AIR_TEMP_50M;
import static psu.ets.psubams.model.SeriesNames.OUT_AIR_TEMP_FOLIAGE;
import static psu.ets.psubams.model.SeriesNames.OUT_BOWEN_RATIO;
import static psu.ets.psubams.model.SeriesNames.OUT_CO2_CONC_CANOPY;
import static psu.ets.psubams.model.SeriesNames.OUT_CO2_FLUX;
import static psu.ets.psubams.model.SeriesNames.OUT_ED1;
import static psu.ets.psubams.model.SeriesNames.OUT_EF1;
import static psu.ets.psubams.model.SeriesNames.OUT_EF2;
import static psu.ets.psubams.model.SeriesNames.OUT_EPIDERMAL_WATER_POT;
import static psu.ets.psubams.model.SeriesNames.OUT_GLOBAL_O3_FLUX;
import static psu.ets.psubams.model.SeriesNames.OUT_GROUND_FLUX;
import static psu.ets.psubams.model.SeriesNames.OUT_GROUND_WATER_POT;
import static psu.ets.psubams.model.SeriesNames.OUT_K_INDEX;
import static psu.ets.psubams.model.SeriesNames.OUT_LATENT_HEAT_FLUX;
import static psu.ets.psubams.model.SeriesNames.OUT_LEAF_WATER_POTENTIAL;
import static psu.ets.psubams.model.SeriesNames.OUT_LFT_INDEX;
import static psu.ets.psubams.model.SeriesNames.OUT_LWDN;
import static psu.ets.psubams.model.SeriesNames.OUT_LWUP;
import static psu.ets.psubams.model.SeriesNames.OUT_NET_RADIATION;
import static psu.ets.psubams.model.SeriesNames.OUT_O3_CONC_CANOPY;
import static psu.ets.psubams.model.SeriesNames.OUT_O3_FLUX_PLANT;
import static psu.ets.psubams.model.SeriesNames.OUT_RADIOMETRIC_TEMP;
import static psu.ets.psubams.model.SeriesNames.OUT_ROOT_ZONE_MOIST_AVAIL;
import static psu.ets.psubams.model.SeriesNames.OUT_SENSIBLE_HEAT_FLUX;
import static psu.ets.psubams.model.SeriesNames.OUT_SHORTWAVE_FLUX;
import static psu.ets.psubams.model.SeriesNames.OUT_SPEC_HUMID_10M;
import static psu.ets.psubams.model.SeriesNames.OUT_SPEC_HUMID_50M;
import static psu.ets.psubams.model.SeriesNames.OUT_SPEC_HUMID_FOLIAGE;
import static psu.ets.psubams.model.SeriesNames.OUT_STOMATAL_RESIST;
import static psu.ets.psubams.model.SeriesNames.OUT_SURF_MOIST_AVAIL;
import static psu.ets.psubams.model.SeriesNames.OUT_SWT_INDEX;
import static psu.ets.psubams.model.SeriesNames.OUT_TIME;
import static psu.ets.psubams.model.SeriesNames.OUT_TOT_TLS;
import static psu.ets.psubams.model.SeriesNames.OUT_VAPOR_PRESS_DEFICIT;
import static psu.ets.psubams.model.SeriesNames.OUT_WATER_USE_EFF;
import static psu.ets.psubams.model.SeriesNames.OUT_WET_BLB_POT_TEMP;
import static psu.ets.psubams.model.SeriesNames.OUT_WIND_10M;
import static psu.ets.psubams.model.SeriesNames.OUT_WIND_50M;
import static psu.ets.psubams.model.SeriesNames.OUT_WIND_IN_FOLIAGE;

public class CompactSimulation implements Constants, SeriesNames {

    public final static int temporal_rows= 24*4+2;
    public static class ModelException extends Exception {

        public ModelException(String s) {
            super(s);
        }
    }

    /**
     * run the model and fill the output arrays
     */
    public static void runModel(Simsphere ps, String csvfile) throws Exception, ModelException, IOException {

        CompactSimulation compactSimulation = new CompactSimulation(ps);        
        String[] output_names = {OUT_TIME, OUT_SHORTWAVE_FLUX, OUT_NET_RADIATION, OUT_SENSIBLE_HEAT_FLUX,
            OUT_LATENT_HEAT_FLUX, OUT_GROUND_FLUX, OUT_AIR_TEMP_50M, OUT_AIR_TEMP_10M,
            OUT_AIR_TEMP_FOLIAGE, OUT_RADIOMETRIC_TEMP, OUT_WIND_50M, OUT_WIND_10M,
            OUT_WIND_IN_FOLIAGE, OUT_SPEC_HUMID_50M, OUT_SPEC_HUMID_10M, OUT_SPEC_HUMID_FOLIAGE,
            OUT_BOWEN_RATIO, OUT_SURF_MOIST_AVAIL, OUT_ROOT_ZONE_MOIST_AVAIL, OUT_STOMATAL_RESIST,
            OUT_VAPOR_PRESS_DEFICIT, OUT_LEAF_WATER_POTENTIAL, OUT_EPIDERMAL_WATER_POT, OUT_GROUND_WATER_POT,
            OUT_CO2_FLUX, OUT_CO2_CONC_CANOPY, OUT_WATER_USE_EFF, OUT_O3_CONC_CANOPY, OUT_GLOBAL_O3_FLUX, OUT_O3_FLUX_PLANT,
            OUT_WET_BLB_POT_TEMP, OUT_LFT_INDEX, OUT_K_INDEX, OUT_TOT_TLS, OUT_SWT_INDEX, OUT_EF1, OUT_EF2, OUT_LWDN, OUT_LWUP, OUT_ED1};

        // use FileWriter constructor that specifies open for appending
        CsvWriter csvOutput = new CsvWriter(new FileWriter(csvfile, false), ',');

        for (String s : output_names) {
            csvOutput.write(s);
        }
        csvOutput.endRecord();

        for (int i = 1; i < compactSimulation.output_data[0].length; i++) {
            double dtime=compactSimulation.output_data[0][i];
            int hour=(int) dtime;
            int minute=(int) Math.round(60*(dtime-hour));
            
            String stime=Integer.toString(hour)+":"+Integer.toString(minute);
            if(hour < 10) {
                stime="0"+stime;
            }
            if(minute == 0) {
                stime=stime+"0";
            }
            csvOutput.write(stime);
            
            for (int j = 1; j <= 34; j++) {
                csvOutput.write(Double.toString(compactSimulation.output_data[j][i]));
            }
            for (int j = 0; j <= 2; j++) {
                csvOutput.write(Double.toString(compactSimulation.data_graphs[j][i]));
            }

            csvOutput.write(Double.toString(compactSimulation.output_data[38][i]));
            csvOutput.write(Double.toString(compactSimulation.output_data[39][i]));
            csvOutput.endRecord();
        }
        csvOutput.close();
    }

    public static void runModelConv(Simsphere ps, String csvfile) throws Exception, ModelException, IOException {               
        int fr_bound= Math.round(1/ps.getOverpass().getSweepFractionalVegetationCover());
        int fi_bound= Math.round(1/ps.getOverpass().getSweepSurfaceMoistureAvail());
        int no_rows=(1+fr_bound)*(1+fi_bound);
        int no_row = 0;
        // Create output arrays
        double[][] output_data = new double[CONVOLUTION_COUNT][no_rows];
        
        int targettime=ps.getOverpass().getRunDate().getHour()*60+ps.getOverpass().getRunDate().getMinute();
        int target_index = targettime / 15;
        for (int i=0;i<=fr_bound;i++) {
            for (int j=0;j<=fi_bound;j++) {
                
                float fr=i*ps.getOverpass().getSweepFractionalVegetationCover();
                if(fr >=1) fr=1;
                if(fr < 0) fr=0;
                float fi=j*ps.getOverpass().getSweepSurfaceMoistureAvail();
                if(fi >=1) fi=1;
                if(fi < 0) fi=0;
                // Call the model     
                
                ps.getVegetation().getVegCommon().setFractVegCover(100*fr);
                ps.getSurface().getHydrologicalParams().setSurfaceMoistureAvail(fi);                
                CompactSimulation compactSimulation = new CompactSimulation(ps); 
                output_data[0][no_row] = 100*fr;      //fractional vegetation cover
                output_data[1][no_row] = fi;   //current surface moisture availability
                output_data[2][no_row] = compactSimulation.output_data[17][target_index] / 1000;   //actual surface moisture availability
                output_data[3][no_row] = compactSimulation.output_data[2][target_index];       //net radiation
                output_data[4][no_row] = compactSimulation.output_data[4][target_index];       //latent heat flux
                output_data[5][no_row] = compactSimulation.output_data[3][target_index];      //sensible heat flux
                output_data[6][no_row] = compactSimulation.output_data[9][target_index]; //radiometric temperature
                no_row++;
            }
        }

        String[] output_names = {"Fr", "Mo", "Mo act", "Rn", "LE",
            "H", "Tir"};

        // use FileWriter constructor that specifies open for appending
        CsvWriter csvOutput = new CsvWriter(new FileWriter(csvfile, false), ',');

        for (String s : output_names) {
            csvOutput.write(s);
        }
        csvOutput.endRecord();

        for (int i = 0; i < output_data[0].length; i++) {
            for (int j = 0; j < output_names.length; j++) {
                csvOutput.write(Double.toString(output_data[j][i]));
            }
            csvOutput.endRecord();
        }
        csvOutput.close();
    }
    
    
    // Create output arrays
    double[][] output_data = new double[OUTPUT_COUNT][temporal_rows];
    double[][][] sounding_data = new double[SOUNDING_COUNT][temporal_rows][SOUNDING_LEVELS];
    double[][] data_graphs = new double[4][temporal_rows];
        
    public CompactSimulation(Simsphere ps) throws Exception, ModelException, IOException {
        SoundingHelper shelper = new SoundingHelper(ps.getSoundingSet());
        

        Plant plant = null;
        String plantSelection = ps.getPlants().getPlantSelection().getID();
        for (Plant temp_plant : ps.getPlants().getPlantSet().getPlant()) {
            if (temp_plant.getID().equals(plantSelection)) {
                plant = temp_plant;
                break;
            }
        }
        Soil soil = null;
        String soilSelection = ps.getSoils().getSoilSelection().getID();
        for (Soil temp_soil : ps.getSoils().getSoilSet().getSoil()) {
            if (temp_soil.getID().equals(soilSelection)) {
                soil = temp_soil;
                break;
            }
        }

        if ((soil == null) || (plant == null)) {
            throw new ModelException("soil or plant null");
        }
        double rmintemp = ps.getVegetation().getVegetationModelSelection().getDeardorff().getBulkStomatalResistance();
        double rcuttemp = ps.getVegetation().getVegetationModelSelection().getDeardorff().getCuticularResistance();

        //Default is Deardorf
        char whichtype = 'D';
        if (ps.getVegetation().getVegetationModelSelection().getModelSelection().value().equals("CarlsonLynn")) {
            whichtype = 'L';
            rmintemp = plant.getMinStomatalResistance();
            rcuttemp = plant.getCuticleResistance();
        } else if (ps.getVegetation().equals("Deardorff")) {
            whichtype = 'D';
        } else {
            throw new ModelException("Ivalid Vegetation Parametrization");
        }

        Model model = new Model();
        // Call the model
        model._psubams(
                ps.getOverpass().getRunDate().getYear(), // iyr
                ps.getOverpass().getRunDate().getMonth(), // imo
                ps.getOverpass().getRunDate().getDay(), //iday
                ps.getOverpass().getRunDate().getTimezone() / 60, // tz
                ps.getOverpass().getLatitude(), // xlat
                ps.getOverpass().getLongitude(), // xlong
                -15 * 60, // strtim
                24 * 3600, // timend
                //			ps.getOverpass().getOutputInterval().getGregorianCalendar().get(GregorianCalendar.MINUTE), // outtt
                15 * 60, // outtt
                ps.getOverpass().getSlope(), // slope
                ps.getOverpass().getAspect(), // aspect
                ps.getSurface().getHydrologicalParams().getSurfaceMoistureAvail(), // f
                ps.getSurface().getHydrologicalParams().getRootZoneMoistureAvail(), // fsub
                ps.getSurface().getHydrologicalParams().getSubstrateMaxVolWaterContent(), //wmax,
                ps.getSurface().getSurfaceParams().getSubstrateClimatologicalMeanTemp(), // btemp
                ps.getSurface().getSurfaceParams().getThermalInertia().getValue(), // tp
                /*update by Yannis Konstas(start)*/
                ps.getSurface().getSurfaceParams().getThermalInertia().isCalculate() ? 'T' : 'F', // dual_ti
                /*update by Yannis Konstas(end)*/
                0.04d, // ti_a,
                0.02d, // ti_b,
                ps.getSurface().getSurfaceParams().getGroundAlbedo().isCalculate() ? 'T' : 'F', // albedo_gflag
                ps.getSurface().getSurfaceParams().getGroundAlbedo().getValue(), // albg
                ps.getSurface().getSurfaceParams().getGroundEmissivity(), // epsi
                ps.getMeteorological().getAtmosphericPrecipitableWater(), // omega
                ps.getMeteorological().getSurfaceRoughness(), // zo
                ps.getMeteorological().getObstacleHeight(), // obst_hgt
                ps.getMeteorological().getCloudCover().isUse() ? 'T' : 'F', // cloud_flag
                ps.getMeteorological().getCloudCover().getValue(), // cld_fract
                ps.getVegetation().getVegCommon().getFractVegCover(), // frveg
                ps.getVegetation().getVegCommon().getLeafAreaIndex(), // xlai
                ps.getVegetation().getVegCommon().getFoliageEmissivity(), // epsf,
                ps.getVegetation().getVegCommon().getFoliageAlbedo().isCalculate() ? 'T' : 'F', // albedo_fflag
                ps.getVegetation().getVegCommon().getFoliageAlbedo().getValue(), // albf
                whichtype, // stmtype
                ps.getVegetation().getVegetationModelSelection().getDeardorff().getRelativeWaterVolume(), // volrel
                rmintemp, //rmin
                rcuttemp, //rcut
                //ps.getVegetation().getBulkStomatalResistance().getValue(), // rmin
                //ps.getVegetation().getCuticularResistance().getValue(), // rcut
                ps.getVegetation().getVegetationModelSelection().getDeardorff().getVolumetricWaterContent(), // wilt
                ps.getVegetation().getVegetationModelSelection().getCarlsonLynn().getVegetationHeight(), // vegheight
                ps.getVegetation().getVegetationModelSelection().getCarlsonLynn().getLeafWidth(), // width
                ps.getVegetation().getVegetationModelSelection().getCarlsonLynn().isNoCapacitance() ? 'T' : 'F', // steady
                ps.getVegetation().getVegCommon().getCI(), // ci
                ps.getVegetation().getVegCommon().getCO(), // co
                ps.getVegetation().getVegCommon().getSurfOzoneConcentration(), // coz_sfc
                ps.getVegetation().getVegCommon().getAirOzoneConcentration(), // coz_air
                soil.getRKS(), // rks
                soil.getCosbyB(), // cosbyb
                soil.getTHMax(), // thmax
                soil.getPsi(), // psis
                plant.getTemperatureRange().getMinTemperature(), // mintemp
                plant.getTemperatureRange().getMaxTemperature(), // maxtemp
                plant.getBeta(), // beta
                plant.getB1(), // b1
                plant.getB2(), // b2
                plant.getCriticalLeafWaterPot(), // psice
                plant.getCriticalSolarParam(), // sc
                plant.getStemResistance(), // zp
                plant.getFractXylemPot(), // getfrhgt
                plant.getFractStemResistance(), // frzp
                plant.getRkocap(), // rkocap
                plant.getRccap(), // rccap
                plant.getRzcap(), // rzcap
                plant.getInitialPlantWaterVol(), // volini
                plant.getZstini(), // zstini
                shelper.len, // nobs_ptq
                shelper.len, // nobs_win 
                ps.getOverpass().getAltitude(), // station_height 
                ps.getSurface().getSurfaceGeostrophicWind().getSGWU(), // ugs 
                ps.getSurface().getSurfaceGeostrophicWind().getSGWV(), // vgs
                shelper.pressure, // ps
                shelper.temperature, // ts HumiditySoundingElement
                shelper.dewpoint, // dep 
                shelper.winddirection, // dd0 
                shelper.windspeed, // ff0 
                shelper.altitude, // zh                 
                data_graphs,
                sounding_data, // sounding_data
                ps.getOverpass().getRunDate().toGregorianCalendar().get(java.util.GregorianCalendar.DAY_OF_YEAR),
                output_data, //output_data
                false //don't override monitor update (one-run diurnal) /*update by Yannis Konstas*/                
        );
    }
    
    
    
    
    

}
