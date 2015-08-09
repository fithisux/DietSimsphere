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
package testxsd;

import java.io.*;
import javax.xml.validation.*;
import javax.xml.transform.stream.*;
import uk.ac.aber.simsphere.*;
import java.io.File;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import psu.ets.psubams.model.*;

/**
 *
 * @author FITHIS
 */
public class TestXSD {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        File xmlFile = new File(args[0]);
        File xsdFile = new File(args[1]);
        String model_csv=args[2];
        String model_conv_csv=args[3];
        
        // 1. Lookup a factory for the W3C XML Schema language
        SchemaFactory factory = SchemaFactory.newInstance("http://www.w3.org/XML/XMLSchema/v1.1");

        try {
            // 2. Compile the schema.
            File schemaLocation = xsdFile;
            Schema schema = factory.newSchema(schemaLocation);

            // 3. Get a validator from the schema.
            Validator validator = schema.newValidator();

            // 4. Parse the document you want to check.
            StreamSource source = new StreamSource(xmlFile);

            // 5. Check the document
            validator.validate(source);
            System.out.println(xmlFile.getName() + " is valid.");
        } catch (org.xml.sax.SAXException ex) {
            System.out.println(xmlFile.getName() + " is not valid because ");
            System.out.println(ex.getMessage());
            System.exit(0);
        } catch (IOException ex) {
            System.out.println(xmlFile.getName() + " is not valid because ");
            System.out.println(ex.getMessage());
            System.exit(0);
        }

        JAXBContext jaxbContext = null;
        Simsphere customer = null;
        try {
            jaxbContext = JAXBContext.newInstance(Simsphere.class);
            Unmarshaller jaxbUnmarshaller = jaxbContext.createUnmarshaller();
            customer = (Simsphere) jaxbUnmarshaller.unmarshal(xmlFile);
        } catch (JAXBException ex) {
            ex.printStackTrace();
            System.out.println(xmlFile.getName() + " cannot unmarshal because ");
            System.out.println(ex.getMessage());
            System.exit(0);
        }

        try {
            Simulation.runModel(customer, model_csv);
            Simulation.runModelConv(customer, model_conv_csv);
        } catch (Simulation.ModelException | IOException ex) {
            ex.printStackTrace();
            System.out.println(xmlFile.getName() + " cannot run computations ");
            System.out.println(ex.getMessage());
            System.exit(0);
        } catch (Exception ex) {
            ex.printStackTrace();
            System.out.println(xmlFile.getName() + " cannot run computations ");
            System.out.println(ex.getMessage());
            System.exit(0);
        }
    }

}
