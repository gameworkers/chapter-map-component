import React, { PureComponent } from "react";
import PropTypes from "prop-types";
import {
  ComposableMap,
  ZoomableGroup,
  Geographies,
  Geography,
  Markers
} from "react-simple-maps";

import mapPoints from "./map_data";

import SvgContentElementWrapperWithDefs from "./SvgContentElementWrapperWithDefs";
import GWUMarker from "./GWUMarker";

let smallWorldDataPromise;
let worldDataPromise;
function getWorldData(pathname) {
  return (worldDataPromise =
    worldDataPromise ||
    fetch(pathname).then(res => res.json())
  );
}
function getSmallWorldData(pathname) {
  return (smallWorldDataPromise =
    smallWorldDataPromise ||
    fetch(pathname).then(res => res.json())
  );
}

const markers = mapPoints
  .filter(point => point.chapter)
  .sort((a, b) => {
    // we want to make sure markers lower on the map are painted in front
    if (a.lat < b.lat) {
      return 1;
    }
    if (a.lat > b.lat) {
      return -1;
    }
    if (a.lng < b.lng) {
      return 1;
    }
    if (a.lng > b.lng) {
      return -1;
    }
    return 0;
  })
  .map(point => {
    return {
      name: point.location,
      coordinates: [point.lng, point.lat]
    };
  });

const countryNameKeys = [
  "ABBREV",
  "CONTINENT",
  "FORMAL_EN",
  "ISO_A2",
  "ISO_A3",
  "NAME",
  "NAME_LONG"
];
function geographyMatchesCountryString(geography, countryString) {
  return countryNameKeys.some(key => {
    return geography.properties[key] === countryString;
  });
}

class ChapterMapComponent extends PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      worldData: null
    };
  }

  componentDidMount() {
    if (
      typeof WORLD_110M_JSON_PATH === 'undefined' ||
      typeof WORLD_50M_JSON_PATH === 'undefined'
    ) {
      throw new Error(`
        WORLD_110M_JSON_PATH and WORLD_50M_JSON_PATH must be defined
        in the global scope.
      `);
    }

    getWorldData(WORLD_50M_JSON_PATH).then(worldData => {
      this.setState({
        worldData
      });
    });
    getSmallWorldData(WORLD_110M_JSON_PATH).then(worldData => {
      this.setState(state => (state.worldData ? null : { worldData }));
    });
  }

  render() {
    const {
      centerLat,
      centerLng,
      width,
      height,
      scale,
      isGeographyIncluded,
      markerScale,
      style,
      className,
      forceGrayscale
    } = this.props;
    const { worldData } = this.state;
    return (
      <div className={className} style={{ ...(style || {}), width, height }}>
        {!worldData && <div>Loading...</div>}
        {worldData && (
          <ComposableMap
            projectionConfig={{ scale: scale }}
            width={width}
            height={height}
            style={{ width: "100%", height: "auto", backgroundColor: '#fff' }}
          >
            <SvgContentElementWrapperWithDefs forceGrayscale={forceGrayscale}>
              <ZoomableGroup center={[centerLng, centerLat]} disablePanning>
                <Geographies geography={this.state.worldData}>
                  {(geographies, projection) =>
                    geographies
                      .filter(isGeographyIncluded)
                      .map((geography, i) => {
                        const hasMatchingPoint = mapPoints.some(point => {
                          return geographyMatchesCountryString(
                            geography,
                            point.country
                          );
                        });
                        const style = {
                          fill: hasMatchingPoint
                            ? "url(#redpattern)"
                            : "url(#hardlyredpattern)",
                          stroke: "#222",
                          strokeWidth: 0.5,
                          outline: "none"
                        };
                        return (
                          <Geography
                            key={i}
                            geography={geography}
                            projection={projection}
                            style={{
                              default: style,
                              hover: style,
                              pressed: style
                            }}
                          />
                        );
                      })
                  }
                </Geographies>
                <Markers>
                  {markers.map(marker => {
                    return (
                      <GWUMarker
                        key={marker.name}
                        marker={marker}
                        scale={markerScale}
                      />
                    );
                  })}
                </Markers>
              </ZoomableGroup>
            </SvgContentElementWrapperWithDefs>
          </ComposableMap>
        )}
      </div>
    );
  }
}

ChapterMapComponent.propTypes = {
  centerLat: PropTypes.number.isRequired,
  centerLng: PropTypes.number.isRequired,
  width: PropTypes.number.isRequired,
  height: PropTypes.number.isRequired,
  scale: PropTypes.number.isRequired,
  isGeographyIncluded: PropTypes.func.isRequired,
  markerScale: PropTypes.number.isRequired,
  forceGrayscale: PropTypes.bool.isRequired
};

ChapterMapComponent.defaultProps = {
  centerLat: 0,
  centerLng: 0,
  width: 980,
  height: 551,
  scale: 205,
  isGeographyIncluded: () => true,
  markerScale: 0.09,
  forceGrayscale: false
};

export default ChapterMapComponent;
