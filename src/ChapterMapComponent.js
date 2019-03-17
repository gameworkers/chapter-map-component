import React, { PureComponent } from "react";
import PropTypes from "prop-types";
import {
  ComposableMap,
  ZoomableGroup,
  Geographies,
  Geography,
  Markers
} from "react-simple-maps";
import withRedux from "next-redux-wrapper";
import { Tooltip, actions } from "redux-tooltip";

React.PropTypes = PropTypes; // for redux-tooltip compatibility

const wrapperStyles = {
  width: "100%",
  maxWidth: 980,
  margin: "0 auto",
  fontFamily: "Roboto, sans-serif"
};

import mapPoints from "./map_data";

import SvgContentElementWrapperWithDefs from "./SvgContentElementWrapperWithDefs";
import GWUMarker from "./GWUMarker";
import { initStore } from "./tooltipStore";

let smallWorldDataPromise;
let worldDataPromise;
function getWorldData(pathname) {
  return (worldDataPromise =
    worldDataPromise || fetch(pathname).then(res => res.json()));
}
function getSmallWorldData(pathname) {
  return (smallWorldDataPromise =
    smallWorldDataPromise || fetch(pathname).then(res => res.json()));
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
      coordinates: [point.lng, point.lat],
      data: point
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
    this.focusedMarker = null;
    this.handleCountryMove = this.handleCountryMove.bind(this);
    this.handleCountryLeave = this.handleCountryLeave.bind(this);
    this.handleMarkerMove = this.handleMarkerMove.bind(this);
    this.handleMarkerLeave = this.handleMarkerLeave.bind(this);
    this.handleMarkerClick = this.handleMarkerClick.bind(this);
    this.handleOutsideMarkerClick = this.handleOutsideMarkerClick.bind(this);
  }

  componentDidMount() {
    if (
      typeof WORLD_110M_JSON_PATH === "undefined" ||
      typeof WORLD_50M_JSON_PATH === "undefined"
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

    window.addEventListener("click", this.handleOutsideMarkerClick);
  }

  componentWillUnmount() {
    window.removeEventListener("click", this.handleOutsideMarkerClick);
  }

  handleCountryMove(geography, evt) {
    if (this.focusedMarker) {
      return;
    }
    this.dispatchTooltip(evt, geography.properties.NAME_LONG);
  }

  handleCountryLeave() {
    if (this.focusedMarker) {
      return;
    }
    this.hideTooltip();
  }

  handleMarkerMove(marker, evt) {
    if (this.focusedMarker) {
      return;
    }
    this.dispatchTooltip(evt, marker.name);
  }

  handleMarkerLeave() {
    if (this.focusedMarker) {
      return;
    }
    this.hideTooltip();
  }

  handleMarkerClick(marker, coords, evt) {
    evt.stopPropagation();
    evt.nativeEvent.stopImmediatePropagation();
    this.focusedMarker = this.focusedMarker === marker ? null : marker;
    if (this.focusedMarker) {
      const {
        data: { location, chapterInfo }
      } = marker;
      let content = location;
      if (chapterInfo) {
        const {
          description,
          applicationLink,
          twitter,
          email,
          website
        } = chapterInfo;
        content = `<h3>${location}</h3>`;
        if (description) {
          content += `<p>${description}</p>`;
        }
        if (twitter || email || website) {
          content += "<p>";
        }
        if (twitter) {
          content += `
            <strong>Twitter:</strong>
            <a href="https://twitter.com/${twitter}">
              @${twitter}
            </a>
            <br />
          `;
        }
        if (email) {
          content += `
            <strong>Email:</strong>
            <a href="mailto:${email}">
              ${email}
            </a>
            <br />
          `;
        }
        if (website) {
          content += `
            <strong>Website:</strong>
            <a href="${
              website.indexOf("http") === 0 ? "" : "http://"
            }${website}">
              ${website}
            </a>
            <br />
          `;
        }
        if (twitter || email || website) {
          content += "</p>";
        }
        if (applicationLink) {
          content += `
            <p>
              <a href="${applicationLink}">
                Apply here!
              </a>
            </p>
          `;
        }
        content = `<div class="${
          this.props.tooltipClassName
        }">${content}</div>`;
      }
      this.dispatchTooltip(evt, content);
    } else {
      this.hideTooltip();
    }
  }

  handleOutsideMarkerClick() {
    this.hideTooltip();
  }

  dispatchTooltip(evt, content) {
    const x = evt.clientX;
    const y = evt.clientY + window.pageYOffset;
    this.props.dispatch(
      actions.show({
        origin: { x, y },
        content
      })
    );
  }

  hideTooltip() {
    this.props.dispatch(actions.hide());
    Promise.resolve().then(() => {
      this.focusedMarker = null;
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
      forceGrayscale,
      zoom,
      enablePanning,
      projection
    } = this.props;
    const { worldData } = this.state;
    return (
      <div className={className} style={{ ...(style || {}), width, height }}>
        {!worldData && <div>Loading...</div>}
        {worldData && (
          <React.Fragment>
            <ComposableMap
              projection={projection}
              projectionConfig={{ scale: scale }}
              width={width}
              height={height}
              style={{ width: "100%", height: "auto", backgroundColor: "#fff" }}
            >
              <SvgContentElementWrapperWithDefs forceGrayscale={forceGrayscale}>
                <ZoomableGroup
                  center={[centerLng, centerLat]}
                  disablePanning={!enablePanning}
                  zoom={zoom}
                >
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
                            strokeWidth: 0.5 / zoom,
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
                              onMouseMove={
                                hasMatchingPoint && this.handleCountryMove
                              }
                              onMouseLeave={
                                hasMatchingPoint && this.handleCountryLeave
                              }
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
                          onClick={this.handleMarkerClick}
                          onMouseMove={this.handleMarkerMove}
                          onMouseLeave={this.handleMarkerLeave}
                        />
                      );
                    })}
                  </Markers>
                </ZoomableGroup>
              </SvgContentElementWrapperWithDefs>
            </ComposableMap>
            <Tooltip />
          </React.Fragment>
        )}
      </div>
    );
  }
}

ChapterMapComponent = withRedux(initStore)(ChapterMapComponent);
ChapterMapComponent.displayName = "ChapterMapComponent";

ChapterMapComponent.propTypes = {
  centerLat: PropTypes.number.isRequired,
  centerLng: PropTypes.number.isRequired,
  width: PropTypes.number.isRequired,
  height: PropTypes.number.isRequired,
  scale: PropTypes.number.isRequired,
  isGeographyIncluded: PropTypes.func.isRequired,
  markerScale: PropTypes.number.isRequired,
  forceGrayscale: PropTypes.bool.isRequired,
  tooltipClassName: PropTypes.string.isRequired,
  zoom: PropTypes.number.isRequired,
  enablePanning: PropTypes.bool.isRequired,
  projection: PropTypes.string
};

ChapterMapComponent.defaultProps = {
  centerLat: 0,
  centerLng: 0,
  width: 980,
  height: 551,
  scale: 205,
  isGeographyIncluded: () => true,
  markerScale: 0.09,
  forceGrayscale: false,
  tooltipClassName: "gwu_chapter_tooltip",
  zoom: 1,
  enablePanning: false,
  projection: "times"
};

export default ChapterMapComponent;
