//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.8-b130911.1802 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2015.03.10 at 10:39:31 AM EET 
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
 *         &lt;element ref="{http://aber.ac.uk/simsphere}SoilSelection"/>
 *         &lt;element ref="{http://aber.ac.uk/simsphere}SoilSet"/>
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
    "soilSelection",
    "soilSet"
})
@XmlRootElement(name = "Soils")
public class Soils {

    @XmlElement(name = "SoilSelection", required = true)
    protected SoilSelection soilSelection;
    @XmlElement(name = "SoilSet", required = true)
    protected SoilSet soilSet;

    /**
     * Gets the value of the soilSelection property.
     * 
     * @return
     *     possible object is
     *     {@link SoilSelection }
     *     
     */
    public SoilSelection getSoilSelection() {
        return soilSelection;
    }

    /**
     * Sets the value of the soilSelection property.
     * 
     * @param value
     *     allowed object is
     *     {@link SoilSelection }
     *     
     */
    public void setSoilSelection(SoilSelection value) {
        this.soilSelection = value;
    }

    /**
     * Gets the value of the soilSet property.
     * 
     * @return
     *     possible object is
     *     {@link SoilSet }
     *     
     */
    public SoilSet getSoilSet() {
        return soilSet;
    }

    /**
     * Sets the value of the soilSet property.
     * 
     * @param value
     *     allowed object is
     *     {@link SoilSet }
     *     
     */
    public void setSoilSet(SoilSet value) {
        this.soilSet = value;
    }

}
