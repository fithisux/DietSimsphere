/*
 * Copyright (C) 2016 orbit
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
package testxsd;
import java.io.FileNotFoundException;
import java.io.IOException;

import com.csvreader.CsvReader;
/**
 *
 * @author orbit
 */
public class CSV2SOUNDINGSET {
    
    public static void main(String[] args){
        
        try {
			
			CsvReader products = new CsvReader(args[0]);
		
			products.readHeaders();

                        int index=0;
                        System.out.println("<SoundingSet>");
			while (products.readRecord())
			{
				String Pressure = products.get("Pressure");
				String Temperature = products.get("Temperature");
				String DewPointDepression = products.get("DewPointDepression");
				String WindSpeed = products.get("WindSpeed");
				String WindDirection = products.get("WindDirection");
				String Height = products.get("Height");
				
                                int a = (int) Math.round(Double.parseDouble(Height) * 3.28084);
                                
                                String AltAboveStation= Integer.toString(a);
                                
                                if(index==0){
                                    AltAboveStation="0";
                                }
                                double b= Math.min(50,Math.abs(Double.parseDouble(DewPointDepression)));
				DewPointDepression= Double.toString(b);
                                System.out.println("<Sounding AltAboveStation=\""+AltAboveStation+"\">");
                                System.out.println("<WindDirection>"+WindDirection+"</WindDirection>");
                                System.out.println("<WindSpeed>"+WindSpeed+"</WindSpeed>");
                                System.out.println("<Pressure>"+Pressure+"</Pressure>");
                                System.out.println("<Temperature>"+Temperature+"</Temperature>");
                                System.out.println("<DewPointDepression>"+DewPointDepression+"</DewPointDepression>");
                                System.out.println("</Sounding>");
                                index++;
                                if(index==40) break;
			}
                        System.out.println("</SoundingSet>");
                        
			products.close();
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
        
    }
    
}
