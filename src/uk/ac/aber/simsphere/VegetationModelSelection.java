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
 *         &lt;element ref="{http://aber.ac.uk/simsphere}Deardorff"/>
 *         &lt;element ref="{http://aber.ac.uk/simsphere}CarlsonLynn"/>
 *       &lt;/sequence>
 *       &lt;attribute name="ModelSelection" use="required" type="{http://aber.ac.uk/simsphere}ModelSelectionType" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "", propOrder = {
    "deardorff",
    "carlsonLynn"
})
@XmlRootElement(name = "VegetationModelSelection")
public class VegetationModelSelection {

    @XmlElement(name = "Deardorff", required = true)
    protected Deardorff deardorff;
    @XmlElement(name = "CarlsonLynn", required = true)
    protected CarlsonLynn carlsonLynn;
    @XmlAttribute(name = "ModelSelection", required = true)
    protected ModelSelectionType modelSelection;

    /**
     * Gets the value of the deardorff property.
     * 
     * @return
     *     possible object is
     *     {@link Deardorff }
     *     
     */
    public Deardorff getDeardorff() {
        return deardorff;
    }

    /**
     * Sets the value of the deardorff property.
     * 
     * @param value
     *     allowed object is
     *     {@link Deardorff }
     *     
     */
    public void setDeardorff(Deardorff value) {
        this.deardorff = value;
    }

    /**
     * Gets the value of the carlsonLynn property.
     * 
     * @return
     *     possible object is
     *     {@link CarlsonLynn }
     *     
     */
    public CarlsonLynn getCarlsonLynn() {
        return carlsonLynn;
    }

    /**
     * Sets the value of the carlsonLynn property.
     * 
     * @param value
     *     allowed object is
     *     {@link CarlsonLynn }
     *     
     */
    public void setCarlsonLynn(CarlsonLynn value) {
        this.carlsonLynn = value;
    }

    /**
     * Gets the value of the modelSelection property.
     * 
     * @return
     *     possible object is
     *     {@link ModelSelectionType }
     *     
     */
    public ModelSelectionType getModelSelection() {
        return modelSelection;
    }

    /**
     * Sets the value of the modelSelection property.
     * 
     * @param value
     *     allowed object is
     *     {@link ModelSelectionType }
     *     
     */
    public void setModelSelection(ModelSelectionType value) {
        this.modelSelection = value;
    }

}
