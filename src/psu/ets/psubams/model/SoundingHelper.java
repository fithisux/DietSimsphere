/*
 * Copyright (C) 2015 Vasileios Anagnostopoulos
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
import uk.ac.aber.simsphere.*;
import java.util.*;

/**
 *
 * @author FITHIS
 */
public class SoundingHelper {
    public int len;
    public double[] pressure;
    public double[] temperature;
    public double[] dewpoint;
    public double[] altitude;
    public double[] windspeed;
    public double[] winddirection;
    
    public SoundingHelper(SoundingSet sset){
        this.len=sset.getSounding().size();
        
        Comparator<Sounding> byAltAboveStation = (e1, e2) -> Integer.compare(
            e1.getAltAboveStation(), e2.getAltAboveStation());

        List<Sounding> sorted=new ArrayList<>();
        sset.getSounding().stream().sorted(byAltAboveStation).forEach(e -> sorted.add(e));
        
        pressure=new double[len+2];
        temperature=new double[len+2];
        dewpoint=new double[len+2];
        altitude=new double[len+1];
        windspeed=new double[len+1];
        winddirection=new double[len+1];
        
        for(int i=0;i<len;i++){
            pressure[i+1]=sorted.get(i).getPressure();
            temperature[i+1]=sorted.get(i).getTemperature();
            dewpoint[i+1]=sorted.get(i).getDewPointDepression();
            altitude[i+1]=( (double) sorted.get(i).getAltAboveStation() ) / 1000 ;
            windspeed[i+1]=sorted.get(i).getWindSpeed();            
            winddirection[i+1]=sorted.get(i).getWindDirection();
        } 
        //pressure[len+1]=250;
        //temperature[len+1]=-40;
        //dewpoint[len+1]=30;
    }
}
