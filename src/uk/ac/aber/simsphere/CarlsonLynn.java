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
 *         &lt;element name="VegetationHeight">
 *           &lt;simpleType>
 *             &lt;restriction base="{http://www.w3.org/2001/XMLSchema}float">
 *               &lt;minInclusive value="0.02"/>
 *               &lt;maxInclusive value="20"/>
 *             &lt;/restriction>
 *           &lt;/simpleType>
 *         &lt;/element>
 *         &lt;element name="LeafWidth">
 *           &lt;simpleType>
 *             &lt;restriction base="{http://www.w3.org/2001/XMLSchema}float">
 *               &lt;minInclusive value="0.01"/>
 *               &lt;maxInclusive value="1"/>
 *             &lt;/restriction>
 *           &lt;/simpleType>
 *         &lt;/element>
 *         &lt;element name="NoCapacitance">
 *           &lt;simpleType>
 *             &lt;restriction base="{http://www.w3.org/2001/XMLSchema}boolean">
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
    "vegetationHeight",
    "leafWidth",
    "noCapacitance"
})
@XmlRootElement(name = "CarlsonLynn")
public class CarlsonLynn {

    @XmlElement(name = "VegetationHeight")
    protected float vegetationHeight;
    @XmlElement(name = "LeafWidth")
    protected float leafWidth;
    @XmlElement(name = "NoCapacitance")
    protected boolean noCapacitance;

    /**
     * Gets the value of the vegetationHeight property.
     * 
     */
    public float getVegetationHeight() {
        return vegetationHeight;
    }

    /**
     * Sets the value of the vegetationHeight property.
     * 
     */
    public void setVegetationHeight(float value) {
        this.vegetationHeight = value;
    }

    /**
     * Gets the value of the leafWidth property.
     * 
     */
    public float getLeafWidth() {
        return leafWidth;
    }

    /**
     * Sets the value of the leafWidth property.
     * 
     */
    public void setLeafWidth(float value) {
        this.leafWidth = value;
    }

    /**
     * Gets the value of the noCapacitance property.
     * 
     */
    public boolean isNoCapacitance() {
        return noCapacitance;
    }

    /**
     * Sets the value of the noCapacitance property.
     * 
     */
    public void setNoCapacitance(boolean value) {
        this.noCapacitance = value;
    }

}
