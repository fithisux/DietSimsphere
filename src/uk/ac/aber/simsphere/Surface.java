//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.8-b130911.1802 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2016.02.10 at 03:01:39 PM EET 
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
 *         &lt;element ref="{http://aber.ac.uk/simsphere}SurfaceParams"/>
 *         &lt;element ref="{http://aber.ac.uk/simsphere}HydrologicalParams"/>
 *         &lt;element ref="{http://aber.ac.uk/simsphere}SurfaceGeostrophicWind"/>
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
    "surfaceParams",
    "hydrologicalParams",
    "surfaceGeostrophicWind"
})
@XmlRootElement(name = "Surface")
public class Surface {

    @XmlElement(name = "SurfaceParams", required = true)
    protected SurfaceParams surfaceParams;
    @XmlElement(name = "HydrologicalParams", required = true)
    protected HydrologicalParams hydrologicalParams;
    @XmlElement(name = "SurfaceGeostrophicWind", required = true)
    protected SurfaceGeostrophicWind surfaceGeostrophicWind;

    /**
     * Gets the value of the surfaceParams property.
     * 
     * @return
     *     possible object is
     *     {@link SurfaceParams }
     *     
     */
    public SurfaceParams getSurfaceParams() {
        return surfaceParams;
    }

    /**
     * Sets the value of the surfaceParams property.
     * 
     * @param value
     *     allowed object is
     *     {@link SurfaceParams }
     *     
     */
    public void setSurfaceParams(SurfaceParams value) {
        this.surfaceParams = value;
    }

    /**
     * Gets the value of the hydrologicalParams property.
     * 
     * @return
     *     possible object is
     *     {@link HydrologicalParams }
     *     
     */
    public HydrologicalParams getHydrologicalParams() {
        return hydrologicalParams;
    }

    /**
     * Sets the value of the hydrologicalParams property.
     * 
     * @param value
     *     allowed object is
     *     {@link HydrologicalParams }
     *     
     */
    public void setHydrologicalParams(HydrologicalParams value) {
        this.hydrologicalParams = value;
    }

    /**
     * Gets the value of the surfaceGeostrophicWind property.
     * 
     * @return
     *     possible object is
     *     {@link SurfaceGeostrophicWind }
     *     
     */
    public SurfaceGeostrophicWind getSurfaceGeostrophicWind() {
        return surfaceGeostrophicWind;
    }

    /**
     * Sets the value of the surfaceGeostrophicWind property.
     * 
     * @param value
     *     allowed object is
     *     {@link SurfaceGeostrophicWind }
     *     
     */
    public void setSurfaceGeostrophicWind(SurfaceGeostrophicWind value) {
        this.surfaceGeostrophicWind = value;
    }

}
