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
 *         &lt;element name="CuticularResistance">
 *           &lt;simpleType>
 *             &lt;restriction base="{http://www.w3.org/2001/XMLSchema}float">
 *               &lt;minInclusive value="50"/>
 *               &lt;maxInclusive value="1000"/>
 *             &lt;/restriction>
 *           &lt;/simpleType>
 *         &lt;/element>
 *         &lt;element name="RelativeWaterVolume">
 *           &lt;simpleType>
 *             &lt;restriction base="{http://www.w3.org/2001/XMLSchema}float">
 *               &lt;minInclusive value="1"/>
 *               &lt;maxInclusive value="100"/>
 *             &lt;/restriction>
 *           &lt;/simpleType>
 *         &lt;/element>
 *         &lt;element name="BulkStomatalResistance">
 *           &lt;simpleType>
 *             &lt;restriction base="{http://www.w3.org/2001/XMLSchema}float">
 *               &lt;minInclusive value="0"/>
 *               &lt;maxInclusive value="250"/>
 *             &lt;/restriction>
 *           &lt;/simpleType>
 *         &lt;/element>
 *         &lt;element name="VolumetricWaterContent">
 *           &lt;simpleType>
 *             &lt;restriction base="{http://www.w3.org/2001/XMLSchema}float">
 *               &lt;minInclusive value="0.05"/>
 *               &lt;maxInclusive value="0.2"/>
 *             &lt;/restriction>
 *           &lt;/simpleType>
 *         &lt;/element>
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
    "cuticularResistance",
    "relativeWaterVolume",
    "bulkStomatalResistance",
    "volumetricWaterContent"
})
@XmlRootElement(name = "Deardorff")
public class Deardorff {

    @XmlElement(name = "CuticularResistance")
    protected float cuticularResistance;
    @XmlElement(name = "RelativeWaterVolume")
    protected float relativeWaterVolume;
    @XmlElement(name = "BulkStomatalResistance")
    protected float bulkStomatalResistance;
    @XmlElement(name = "VolumetricWaterContent")
    protected float volumetricWaterContent;

    /**
     * Gets the value of the cuticularResistance property.
     * 
     */
    public float getCuticularResistance() {
        return cuticularResistance;
    }

    /**
     * Sets the value of the cuticularResistance property.
     * 
     */
    public void setCuticularResistance(float value) {
        this.cuticularResistance = value;
    }

    /**
     * Gets the value of the relativeWaterVolume property.
     * 
     */
    public float getRelativeWaterVolume() {
        return relativeWaterVolume;
    }

    /**
     * Sets the value of the relativeWaterVolume property.
     * 
     */
    public void setRelativeWaterVolume(float value) {
        this.relativeWaterVolume = value;
    }

    /**
     * Gets the value of the bulkStomatalResistance property.
     * 
     */
    public float getBulkStomatalResistance() {
        return bulkStomatalResistance;
    }

    /**
     * Sets the value of the bulkStomatalResistance property.
     * 
     */
    public void setBulkStomatalResistance(float value) {
        this.bulkStomatalResistance = value;
    }

    /**
     * Gets the value of the volumetricWaterContent property.
     * 
     */
    public float getVolumetricWaterContent() {
        return volumetricWaterContent;
    }

    /**
     * Sets the value of the volumetricWaterContent property.
     * 
     */
    public void setVolumetricWaterContent(float value) {
        this.volumetricWaterContent = value;
    }

}
