//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.8-b130911.1802 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2015.03.10 at 10:39:31 AM EET 
//


package uk.ac.aber.simsphere;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;
import javax.xml.bind.annotation.XmlValue;


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
 *         &lt;element name="SubstrateClimatologicalMeanTemp">
 *           &lt;simpleType>
 *             &lt;restriction base="{http://www.w3.org/2001/XMLSchema}float">
 *               &lt;minInclusive value="-50"/>
 *               &lt;maxInclusive value="120"/>
 *             &lt;/restriction>
 *           &lt;/simpleType>
 *         &lt;/element>
 *         &lt;element name="ThermalInertia">
 *           &lt;complexType>
 *             &lt;simpleContent>
 *               &lt;extension base="&lt;http://aber.ac.uk/simsphere>ThermalInertiaType">
 *                 &lt;attribute name="calculate" use="required" type="{http://www.w3.org/2001/XMLSchema}boolean" />
 *               &lt;/extension>
 *             &lt;/simpleContent>
 *           &lt;/complexType>
 *         &lt;/element>
 *         &lt;element name="GroundAlbedo">
 *           &lt;complexType>
 *             &lt;simpleContent>
 *               &lt;extension base="&lt;http://aber.ac.uk/simsphere>GroundAlbedoType">
 *                 &lt;attribute name="calculate" use="required" type="{http://www.w3.org/2001/XMLSchema}boolean" />
 *               &lt;/extension>
 *             &lt;/simpleContent>
 *           &lt;/complexType>
 *         &lt;/element>
 *         &lt;element name="GroundEmissivity">
 *           &lt;simpleType>
 *             &lt;restriction base="{http://www.w3.org/2001/XMLSchema}float">
 *               &lt;minInclusive value="0.95"/>
 *               &lt;maxInclusive value="1"/>
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
    "substrateClimatologicalMeanTemp",
    "thermalInertia",
    "groundAlbedo",
    "groundEmissivity"
})
@XmlRootElement(name = "SurfaceParams")
public class SurfaceParams {

    @XmlElement(name = "SubstrateClimatologicalMeanTemp")
    protected float substrateClimatologicalMeanTemp;
    @XmlElement(name = "ThermalInertia", required = true)
    protected SurfaceParams.ThermalInertia thermalInertia;
    @XmlElement(name = "GroundAlbedo", required = true)
    protected SurfaceParams.GroundAlbedo groundAlbedo;
    @XmlElement(name = "GroundEmissivity")
    protected float groundEmissivity;

    /**
     * Gets the value of the substrateClimatologicalMeanTemp property.
     * 
     */
    public float getSubstrateClimatologicalMeanTemp() {
        return substrateClimatologicalMeanTemp;
    }

    /**
     * Sets the value of the substrateClimatologicalMeanTemp property.
     * 
     */
    public void setSubstrateClimatologicalMeanTemp(float value) {
        this.substrateClimatologicalMeanTemp = value;
    }

    /**
     * Gets the value of the thermalInertia property.
     * 
     * @return
     *     possible object is
     *     {@link SurfaceParams.ThermalInertia }
     *     
     */
    public SurfaceParams.ThermalInertia getThermalInertia() {
        return thermalInertia;
    }

    /**
     * Sets the value of the thermalInertia property.
     * 
     * @param value
     *     allowed object is
     *     {@link SurfaceParams.ThermalInertia }
     *     
     */
    public void setThermalInertia(SurfaceParams.ThermalInertia value) {
        this.thermalInertia = value;
    }

    /**
     * Gets the value of the groundAlbedo property.
     * 
     * @return
     *     possible object is
     *     {@link SurfaceParams.GroundAlbedo }
     *     
     */
    public SurfaceParams.GroundAlbedo getGroundAlbedo() {
        return groundAlbedo;
    }

    /**
     * Sets the value of the groundAlbedo property.
     * 
     * @param value
     *     allowed object is
     *     {@link SurfaceParams.GroundAlbedo }
     *     
     */
    public void setGroundAlbedo(SurfaceParams.GroundAlbedo value) {
        this.groundAlbedo = value;
    }

    /**
     * Gets the value of the groundEmissivity property.
     * 
     */
    public float getGroundEmissivity() {
        return groundEmissivity;
    }

    /**
     * Sets the value of the groundEmissivity property.
     * 
     */
    public void setGroundEmissivity(float value) {
        this.groundEmissivity = value;
    }


    /**
     * <p>Java class for anonymous complex type.
     * 
     * <p>The following schema fragment specifies the expected content contained within this class.
     * 
     * <pre>
     * &lt;complexType>
     *   &lt;simpleContent>
     *     &lt;extension base="&lt;http://aber.ac.uk/simsphere>GroundAlbedoType">
     *       &lt;attribute name="calculate" use="required" type="{http://www.w3.org/2001/XMLSchema}boolean" />
     *     &lt;/extension>
     *   &lt;/simpleContent>
     * &lt;/complexType>
     * </pre>
     * 
     * 
     */
    @XmlAccessorType(XmlAccessType.FIELD)
    @XmlType(name = "", propOrder = {
        "value"
    })
    public static class GroundAlbedo {

        @XmlValue
        protected float value;
        @XmlAttribute(name = "calculate", required = true)
        protected boolean calculate;

        /**
         * Gets the value of the value property.
         * 
         */
        public float getValue() {
            return value;
        }

        /**
         * Sets the value of the value property.
         * 
         */
        public void setValue(float value) {
            this.value = value;
        }

        /**
         * Gets the value of the calculate property.
         * 
         */
        public boolean isCalculate() {
            return calculate;
        }

        /**
         * Sets the value of the calculate property.
         * 
         */
        public void setCalculate(boolean value) {
            this.calculate = value;
        }

    }


    /**
     * <p>Java class for anonymous complex type.
     * 
     * <p>The following schema fragment specifies the expected content contained within this class.
     * 
     * <pre>
     * &lt;complexType>
     *   &lt;simpleContent>
     *     &lt;extension base="&lt;http://aber.ac.uk/simsphere>ThermalInertiaType">
     *       &lt;attribute name="calculate" use="required" type="{http://www.w3.org/2001/XMLSchema}boolean" />
     *     &lt;/extension>
     *   &lt;/simpleContent>
     * &lt;/complexType>
     * </pre>
     * 
     * 
     */
    @XmlAccessorType(XmlAccessType.FIELD)
    @XmlType(name = "", propOrder = {
        "value"
    })
    public static class ThermalInertia {

        @XmlValue
        protected float value;
        @XmlAttribute(name = "calculate", required = true)
        protected boolean calculate;

        /**
         * Gets the value of the value property.
         * 
         */
        public float getValue() {
            return value;
        }

        /**
         * Sets the value of the value property.
         * 
         */
        public void setValue(float value) {
            this.value = value;
        }

        /**
         * Gets the value of the calculate property.
         * 
         */
        public boolean isCalculate() {
            return calculate;
        }

        /**
         * Sets the value of the calculate property.
         * 
         */
        public void setCalculate(boolean value) {
            this.calculate = value;
        }

    }

}