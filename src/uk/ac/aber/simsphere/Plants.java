//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.8-b130911.1802 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2017.11.23 at 06:20:10 PM EET 
//


package uk.ac.aber.simsphere;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java class for anonymous complex type.
 * 
 * <p>The following schema fragment specifies the expected content contained within this class.
 * 
 * <pre>
 * &lt;complexType>
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;sequence>
 *         &lt;element ref="{http://aber.ac.uk/simsphere}PlantSelection"/>
 *         &lt;element ref="{http://aber.ac.uk/simsphere}PlantSet"/>
 *       &lt;/sequence>
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "", propOrder = {
    "plantSelection",
    "plantSet"
})
@XmlRootElement(name = "Plants")
public class Plants {

    @XmlElement(name = "PlantSelection", required = true)
    protected PlantSelection plantSelection;
    @XmlElement(name = "PlantSet", required = true)
    protected PlantSet plantSet;

    /**
     * Gets the value of the plantSelection property.
     * 
     * @return
     *     possible object is
     *     {@link PlantSelection }
     *     
     */
    public PlantSelection getPlantSelection() {
        return plantSelection;
    }

    /**
     * Sets the value of the plantSelection property.
     * 
     * @param value
     *     allowed object is
     *     {@link PlantSelection }
     *     
     */
    public void setPlantSelection(PlantSelection value) {
        this.plantSelection = value;
    }

    /**
     * Gets the value of the plantSet property.
     * 
     * @return
     *     possible object is
     *     {@link PlantSet }
     *     
     */
    public PlantSet getPlantSet() {
        return plantSet;
    }

    /**
     * Sets the value of the plantSet property.
     * 
     * @param value
     *     allowed object is
     *     {@link PlantSet }
     *     
     */
    public void setPlantSet(PlantSet value) {
        this.plantSet = value;
    }

}
