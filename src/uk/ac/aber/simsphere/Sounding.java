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
 *         &lt;element name="WindDirection">
 *           &lt;simpleType>
 *             &lt;restriction base="{http://www.w3.org/2001/XMLSchema}float">
 *               &lt;minInclusive value="0"/>
 *               &lt;maxInclusive value="360"/>
 *             &lt;/restriction>
 *           &lt;/simpleType>
 *         &lt;/element>
 *         &lt;element name="WindSpeed">
 *           &lt;simpleType>
 *             &lt;restriction base="{http://www.w3.org/2001/XMLSchema}float">
 *               &lt;minInclusive value="0"/>
 *               &lt;maxInclusive value="100"/>
 *             &lt;/restriction>
 *           &lt;/simpleType>
 *         &lt;/element>
 *         &lt;element name="Pressure">
 *           &lt;simpleType>
 *             &lt;restriction base="{http://www.w3.org/2001/XMLSchema}float">
 *               &lt;minInclusive value="0"/>
 *               &lt;maxInclusive value="1100"/>
 *             &lt;/restriction>
 *           &lt;/simpleType>
 *         &lt;/element>
 *         &lt;element name="Temperature">
 *           &lt;simpleType>
 *             &lt;restriction base="{http://www.w3.org/2001/XMLSchema}float">
 *               &lt;minInclusive value="-100"/>
 *               &lt;maxInclusive value="100"/>
 *             &lt;/restriction>
 *           &lt;/simpleType>
 *         &lt;/element>
 *         &lt;element name="DewPointDepression">
 *           &lt;simpleType>
 *             &lt;restriction base="{http://www.w3.org/2001/XMLSchema}float">
 *               &lt;minInclusive value="0"/>
 *               &lt;maxInclusive value="50"/>
 *             &lt;/restriction>
 *           &lt;/simpleType>
 *         &lt;/element>
 *       &lt;/sequence>
 *       &lt;attribute name="AltAboveStation" type="{http://aber.ac.uk/simsphere}AltAboveStationType" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "", propOrder = {
    "windDirection",
    "windSpeed",
    "pressure",
    "temperature",
    "dewPointDepression"
})
@XmlRootElement(name = "Sounding")
public class Sounding {

    @XmlElement(name = "WindDirection")
    protected float windDirection;
    @XmlElement(name = "WindSpeed")
    protected float windSpeed;
    @XmlElement(name = "Pressure")
    protected float pressure;
    @XmlElement(name = "Temperature")
    protected float temperature;
    @XmlElement(name = "DewPointDepression")
    protected float dewPointDepression;
    @XmlAttribute(name = "AltAboveStation")
    protected Integer altAboveStation;

    /**
     * Gets the value of the windDirection property.
     * 
     */
    public float getWindDirection() {
        return windDirection;
    }

    /**
     * Sets the value of the windDirection property.
     * 
     */
    public void setWindDirection(float value) {
        this.windDirection = value;
    }

    /**
     * Gets the value of the windSpeed property.
     * 
     */
    public float getWindSpeed() {
        return windSpeed;
    }

    /**
     * Sets the value of the windSpeed property.
     * 
     */
    public void setWindSpeed(float value) {
        this.windSpeed = value;
    }

    /**
     * Gets the value of the pressure property.
     * 
     */
    public float getPressure() {
        return pressure;
    }

    /**
     * Sets the value of the pressure property.
     * 
     */
    public void setPressure(float value) {
        this.pressure = value;
    }

    /**
     * Gets the value of the temperature property.
     * 
     */
    public float getTemperature() {
        return temperature;
    }

    /**
     * Sets the value of the temperature property.
     * 
     */
    public void setTemperature(float value) {
        this.temperature = value;
    }

    /**
     * Gets the value of the dewPointDepression property.
     * 
     */
    public float getDewPointDepression() {
        return dewPointDepression;
    }

    /**
     * Sets the value of the dewPointDepression property.
     * 
     */
    public void setDewPointDepression(float value) {
        this.dewPointDepression = value;
    }

    /**
     * Gets the value of the altAboveStation property.
     * 
     * @return
     *     possible object is
     *     {@link Integer }
     *     
     */
    public Integer getAltAboveStation() {
        return altAboveStation;
    }

    /**
     * Sets the value of the altAboveStation property.
     * 
     * @param value
     *     allowed object is
     *     {@link Integer }
     *     
     */
    public void setAltAboveStation(Integer value) {
        this.altAboveStation = value;
    }

}
