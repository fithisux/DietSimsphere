//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.8-b130911.1802 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2016.02.10 at 03:01:39 PM EET 
//


package uk.ac.aber.simsphere;

import javax.xml.bind.annotation.XmlEnum;
import javax.xml.bind.annotation.XmlEnumValue;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java class for ModelSelectionType.
 * 
 * <p>The following schema fragment specifies the expected content contained within this class.
 * <p>
 * <pre>
 * &lt;simpleType name="ModelSelectionType">
 *   &lt;restriction base="{http://www.w3.org/2001/XMLSchema}string">
 *     &lt;enumeration value="Deardorff"/>
 *     &lt;enumeration value="CarlsonLynn"/>
 *   &lt;/restriction>
 * &lt;/simpleType>
 * </pre>
 * 
 */
@XmlType(name = "ModelSelectionType")
@XmlEnum
public enum ModelSelectionType {

    @XmlEnumValue("Deardorff")
    DEARDORFF("Deardorff"),
    @XmlEnumValue("CarlsonLynn")
    CARLSON_LYNN("CarlsonLynn");
    private final String value;

    ModelSelectionType(String v) {
        value = v;
    }

    public String value() {
        return value;
    }

    public static ModelSelectionType fromValue(String v) {
        for (ModelSelectionType c: ModelSelectionType.values()) {
            if (c.value.equals(v)) {
                return c;
            }
        }
        throw new IllegalArgumentException(v);
    }

}
