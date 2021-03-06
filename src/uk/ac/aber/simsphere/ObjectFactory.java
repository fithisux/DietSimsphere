//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.8-b130911.1802 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2017.11.23 at 06:20:10 PM EET 
//


package uk.ac.aber.simsphere;

import javax.xml.bind.JAXBElement;
import javax.xml.bind.annotation.XmlElementDecl;
import javax.xml.bind.annotation.XmlRegistry;
import javax.xml.namespace.QName;


/**
 * This object contains factory methods for each 
 * Java content interface and Java element interface 
 * generated in the uk.ac.aber.simsphere package. 
 * <p>An ObjectFactory allows you to programatically 
 * construct new instances of the Java representation 
 * for XML content. The Java representation of XML 
 * content can consist of schema derived interfaces 
 * and classes representing the binding of schema 
 * type definitions, element declarations and model 
 * groups.  Factory methods for each of these are 
 * provided in this class.
 * 
 */
@XmlRegistry
public class ObjectFactory {

    private final static QName _ObstacleHeight_QNAME = new QName("http://aber.ac.uk/simsphere", "ObstacleHeight");
    private final static QName _SurfaceRoughness_QNAME = new QName("http://aber.ac.uk/simsphere", "SurfaceRoughness");
    private final static QName _AtmosphericPrecipitableWater_QNAME = new QName("http://aber.ac.uk/simsphere", "AtmosphericPrecipitableWater");

    /**
     * Create a new ObjectFactory that can be used to create new instances of schema derived classes for package: uk.ac.aber.simsphere
     * 
     */
    public ObjectFactory() {
    }

    /**
     * Create an instance of {@link SurfaceParams }
     * 
     */
    public SurfaceParams createSurfaceParams() {
        return new SurfaceParams();
    }

    /**
     * Create an instance of {@link VegCommon }
     * 
     */
    public VegCommon createVegCommon() {
        return new VegCommon();
    }

    /**
     * Create an instance of {@link HydrologicalParams }
     * 
     */
    public HydrologicalParams createHydrologicalParams() {
        return new HydrologicalParams();
    }

    /**
     * Create an instance of {@link SoilSet }
     * 
     */
    public SoilSet createSoilSet() {
        return new SoilSet();
    }

    /**
     * Create an instance of {@link Soil }
     * 
     */
    public Soil createSoil() {
        return new Soil();
    }

    /**
     * Create an instance of {@link Deardorff }
     * 
     */
    public Deardorff createDeardorff() {
        return new Deardorff();
    }

    /**
     * Create an instance of {@link VegetationModelSelection }
     * 
     */
    public VegetationModelSelection createVegetationModelSelection() {
        return new VegetationModelSelection();
    }

    /**
     * Create an instance of {@link CarlsonLynn }
     * 
     */
    public CarlsonLynn createCarlsonLynn() {
        return new CarlsonLynn();
    }

    /**
     * Create an instance of {@link Surface }
     * 
     */
    public Surface createSurface() {
        return new Surface();
    }

    /**
     * Create an instance of {@link SurfaceParams.ThermalInertia }
     * 
     */
    public SurfaceParams.ThermalInertia createSurfaceParamsThermalInertia() {
        return new SurfaceParams.ThermalInertia();
    }

    /**
     * Create an instance of {@link SurfaceParams.GroundAlbedo }
     * 
     */
    public SurfaceParams.GroundAlbedo createSurfaceParamsGroundAlbedo() {
        return new SurfaceParams.GroundAlbedo();
    }

    /**
     * Create an instance of {@link SurfaceGeostrophicWind }
     * 
     */
    public SurfaceGeostrophicWind createSurfaceGeostrophicWind() {
        return new SurfaceGeostrophicWind();
    }

    /**
     * Create an instance of {@link Plant }
     * 
     */
    public Plant createPlant() {
        return new Plant();
    }

    /**
     * Create an instance of {@link PlantType }
     * 
     */
    public PlantType createPlantType() {
        return new PlantType();
    }

    /**
     * Create an instance of {@link TemperatureRange }
     * 
     */
    public TemperatureRange createTemperatureRange() {
        return new TemperatureRange();
    }

    /**
     * Create an instance of {@link SoilSelection }
     * 
     */
    public SoilSelection createSoilSelection() {
        return new SoilSelection();
    }

    /**
     * Create an instance of {@link CloudCover }
     * 
     */
    public CloudCover createCloudCover() {
        return new CloudCover();
    }

    /**
     * Create an instance of {@link VegCommon.FoliageAlbedo }
     * 
     */
    public VegCommon.FoliageAlbedo createVegCommonFoliageAlbedo() {
        return new VegCommon.FoliageAlbedo();
    }

    /**
     * Create an instance of {@link PlantSelection }
     * 
     */
    public PlantSelection createPlantSelection() {
        return new PlantSelection();
    }

    /**
     * Create an instance of {@link PlantSet }
     * 
     */
    public PlantSet createPlantSet() {
        return new PlantSet();
    }

    /**
     * Create an instance of {@link Sounding }
     * 
     */
    public Sounding createSounding() {
        return new Sounding();
    }

    /**
     * Create an instance of {@link Plants }
     * 
     */
    public Plants createPlants() {
        return new Plants();
    }

    /**
     * Create an instance of {@link SoundingSet }
     * 
     */
    public SoundingSet createSoundingSet() {
        return new SoundingSet();
    }

    /**
     * Create an instance of {@link Simsphere }
     * 
     */
    public Simsphere createSimsphere() {
        return new Simsphere();
    }

    /**
     * Create an instance of {@link Overpass }
     * 
     */
    public Overpass createOverpass() {
        return new Overpass();
    }

    /**
     * Create an instance of {@link Vegetation }
     * 
     */
    public Vegetation createVegetation() {
        return new Vegetation();
    }

    /**
     * Create an instance of {@link Meteorological }
     * 
     */
    public Meteorological createMeteorological() {
        return new Meteorological();
    }

    /**
     * Create an instance of {@link Soils }
     * 
     */
    public Soils createSoils() {
        return new Soils();
    }

    /**
     * Create an instance of {@link JAXBElement }{@code <}{@link Float }{@code >}}
     * 
     */
    @XmlElementDecl(namespace = "http://aber.ac.uk/simsphere", name = "ObstacleHeight")
    public JAXBElement<Float> createObstacleHeight(Float value) {
        return new JAXBElement<Float>(_ObstacleHeight_QNAME, Float.class, null, value);
    }

    /**
     * Create an instance of {@link JAXBElement }{@code <}{@link Float }{@code >}}
     * 
     */
    @XmlElementDecl(namespace = "http://aber.ac.uk/simsphere", name = "SurfaceRoughness")
    public JAXBElement<Float> createSurfaceRoughness(Float value) {
        return new JAXBElement<Float>(_SurfaceRoughness_QNAME, Float.class, null, value);
    }

    /**
     * Create an instance of {@link JAXBElement }{@code <}{@link Float }{@code >}}
     * 
     */
    @XmlElementDecl(namespace = "http://aber.ac.uk/simsphere", name = "AtmosphericPrecipitableWater")
    public JAXBElement<Float> createAtmosphericPrecipitableWater(Float value) {
        return new JAXBElement<Float>(_AtmosphericPrecipitableWater_QNAME, Float.class, null, value);
    }

}
