//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.8-b130911.1802 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2017.11.23 at 06:20:10 PM EET 
//


package uk.ac.aber.simsphere;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
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
 *         &lt;element name="SGW_U">
 *           &lt;simpleType>
 *             &lt;restriction base="{http://www.w3.org/2001/XMLSchema}float">
 *               &lt;minInclusive value="-30"/>
 *               &lt;maxInclusive value="30"/>
 *             &lt;/restriction>
 *           &lt;/simpleType>
 *         &lt;/element>
 *         &lt;element name="SGW_V">
 *           &lt;simpleType>
 *             &lt;restriction base="{http://www.w3.org/2001/XMLSchema}float">
 *               &lt;minInclusive value="-30"/>
 *               &lt;maxInclusive value="30"/>
 *             &lt;/restriction>
 *           &lt;/simpleType>
 *         &lt;/element>
 *       &lt;/sequence>
 *       &lt;attribute name="UseDefault" use="required" type="{http://www.w3.org/2001/XMLSchema}boolean" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "", propOrder = {
    "sgwu",
    "sgwv"
})
@XmlRootElement(name = "SurfaceGeostrophicWind")
public class SurfaceGeostrophicWind {

    @XmlElement(name = "SGW_U")
    protected float sgwu;
    @XmlElement(name = "SGW_V")
    protected float sgwv;
    @XmlAttribute(name = "UseDefault", required = true)
    protected boolean useDefault;

    /**
     * Gets the value of the sgwu property.
     * 
     */
    public float getSGWU() {
        return sgwu;
    }

    /**
     * Sets the value of the sgwu property.
     * 
     */
    public void setSGWU(float value) {
        this.sgwu = value;
    }

    /**
     * Gets the value of the sgwv property.
     * 
     */
    public float getSGWV() {
        return sgwv;
    }

    /**
     * Sets the value of the sgwv property.
     * 
     */
    public void setSGWV(float value) {
        this.sgwv = value;
    }

    /**
     * Gets the value of the useDefault property.
     * 
     */
    public boolean isUseDefault() {
        return useDefault;
    }

    /**
     * Sets the value of the useDefault property.
     * 
     */
    public void setUseDefault(boolean value) {
        this.useDefault = value;
    }

}
