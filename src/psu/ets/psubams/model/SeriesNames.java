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
package psu.ets.psubams.model;

/**
	constants for names of the different model output series
*/
public interface SeriesNames {

    public static final String OUT_TIME = "Time";
    public static final String OUT_SHORTWAVE_FLUX = "Shortwave Flux (Wm-2)";
    public static final String OUT_NET_RADIATION = "Net Radiation (Wm-2)";
    public static final String OUT_SENSIBLE_HEAT_FLUX = "Sensible Heat Flux (Wm-2)";
    public static final String OUT_LATENT_HEAT_FLUX = "Latent Heat Flux (Wm-2)";
    public static final String OUT_GROUND_FLUX = "Ground Flux (Wm-2)";
    public static final String OUT_AIR_TEMP_50M = "Air Temperature at 50m (C)";
    public static final String OUT_AIR_TEMP_10M = "Air Temperature at 1.3m (C)";
    public static final String OUT_AIR_TEMP_FOLIAGE = "Air Temperature in Foliage (C)";
    public static final String OUT_RADIOMETRIC_TEMP = "Radiometric Temperature (C)";
    public static final String OUT_WIND_50M = "Wind at 50 m (kts)";
    public static final String OUT_WIND_10M = "Wind at 10 m (kts)";
    public static final String OUT_WIND_IN_FOLIAGE = "Wind in Foliage (kts)";
    public static final String OUT_SPEC_HUMID_50M = "Specific Humidity at 50m (gKg-1)";
    public static final String OUT_SPEC_HUMID_10M = "specific Humidity at 1.3m (gKg-1)";
    public static final String OUT_SPEC_HUMID_FOLIAGE = "Specific Humidity of Foliage (gKg-1)";
    public static final String OUT_BOWEN_RATIO = "Bowen Ratio ( x 1000)";
    public static final String OUT_SURF_MOIST_AVAIL = "Surface Moisture Availability ( x 1000)";
    public static final String OUT_ROOT_ZONE_MOIST_AVAIL = "Root Zone Moisture Availability ( x 1000)";
    public static final String OUT_STOMATAL_RESIST = "Stomatal Resistance (sm-1)";
    public static final String OUT_VAPOR_PRESS_DEFICIT = "Vapor Pressure Deficit (mbar)";
    public static final String OUT_LEAF_WATER_POTENTIAL = "Leaf Water Potential (bars)";
    public static final String OUT_EPIDERMAL_WATER_POT = "Epidermal Water Potential (bars)";
    public static final String OUT_GROUND_WATER_POT = "Ground Water Potential (mbars)";
    public static final String OUT_CO2_FLUX = "CO2 Flux (micromolesm-2s-1)";
    public static final String OUT_CO2_CONC_CANOPY = "CO2 Concentration in Canopy (ppmv)";
    public static final String OUT_WATER_USE_EFF = "Water Use Efficiency ( x 1000)";
    public static final String OUT_O3_CONC_CANOPY = "O3 Concentration in Canopy (ppmv ( x 1000))";
    public static final String OUT_GLOBAL_O3_FLUX = "Global O3 Flux (ugm-2s-1)";
    public static final String OUT_O3_FLUX_PLANT = "O3 Flux in Plant (ugm-2s-1)";
    
    // New parameters introduced in New Indices
    public static final String OUT_WET_BLB_POT_TEMP = "Wet Bulb Potential Temperature (K)";
    public static final String OUT_LFT_INDEX = "Lifted Index";
    public static final String OUT_K_INDEX = "K Index";
    public static final String OUT_TOT_TLS = "Total Totals";
    public static final String OUT_SWT_INDEX = "Sweat Index";

    
    /* new parameters */
    //public static final String OUT_HEAT_FLUX_LEAVES = "Heat Flux From Leaves";
    
    public static final String OUT_SOUNDING_HEIGHT = "Height";
    public static final String OUT_SOUNDING_PRESSURE = "Pressure";
    public static final String OUT_SOUNDING_TEMPERATURE = "Temperature";
    public static final String OUT_SOUNDING_HUMIDITY = "Humidity";
    public static final String OUT_SOUNDING_WIND_DIRECTION = "Wind direction";
    public static final String OUT_SOUNDING_WIND_SPEED = "Wind speed";
    public static final String OUT_SOUNDING_POT_TEMP = "Potential temperature";

        /*Update by Yannis Konstas(start)*/
    public static final String OUT_FRACT_VEG_COVER = "Fractional vegetation cover (%)";
    public static final String OUT_EF1 = "LE _ (LE + H)";
    public static final String OUT_EF2 = "H _ (LE + H)";
    public static final String OUT_ED1 = "Daily Heat Flux";
    public static final String OUT_LWUP = "Longwave upwelling radiation";
    public static final String OUT_LWDN = "Longwave downwelling radiation";
    public static final String OUT_SWAVE_UP = "Shortwave Upwelling";

}