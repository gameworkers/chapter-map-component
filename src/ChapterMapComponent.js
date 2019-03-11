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

import GWUMarker from "./GWUMarker";

let smallWorldDataPromise;
let worldDataPromise;
function getWorldData() {
  return (worldDataPromise =
    worldDataPromise || fetch(location.pathname + "/dist/world-50m.json").then(res => res.json()));
}
function getSmallWorldData() {
  return (smallWorldDataPromise =
    smallWorldDataPromise || fetch(location.pathname + "/distworld-110m.json").then(res => res.json()));
}

const markers = mapPoints
  .filter(point => point.chapter)
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
    getWorldData().then(worldData => {
      this.setState({
        worldData
      });
    });
    getSmallWorldData().then(worldData => {
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
      className
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
            style={{ width: "100%", height: "auto" }}
          >
            <ZoomableGroup center={[centerLng, centerLat]} disablePanning>
              <Geographies geography={this.state.worldData}>
                {(geographies, projection) =>
                  geographies
                    .filter(isGeographyIncluded)
                    .map((geography, i) => {
                      const isHighlighted = mapPoints.some(point => {
                        return geographyMatchesCountryString(
                          geography,
                          point.country
                        );
                      });
                      return (
                        <Geography
                          key={i}
                          geography={geography}
                          projection={projection}
                          style={{
                            default: {
                              fill: isHighlighted ? "#ccc" : "#efefef",
                              stroke: "#333",
                              strokeWidth: 0.3,
                              outline: "none"
                            },
                            hover: {
                              fill: "#607D8B",
                              stroke: "#607D8B",
                              strokeWidth: 0.75,
                              outline: "none"
                            },
                            pressed: {
                              fill: "#FF5722",
                              stroke: "#607D8B",
                              strokeWidth: 0.75,
                              outline: "none"
                            }
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
  markerScale: PropTypes.number.isRequired
};

ChapterMapComponent.defaultProps = {
  centerLat: 0,
  centerLng: 0,
  width: 980,
  height: 551,
  scale: 205,
  isGeographyIncluded: () => true,
  markerScale: 0.1
};

export default ChapterMapComponent;
