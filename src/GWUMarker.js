import React, { PureComponent } from "react";
import PropTypes from "prop-types";
import { Marker } from "@gameworkers/react-simple-maps";

class GWUMarker extends PureComponent {
  render() {
    const { marker, scale, ...rest } = this.props;
    return (
      <Marker key={marker.name} marker={marker} {...rest}>
        <g
          style={{ cursor: "pointer" }}
          transform={`scale(${scale}), translate(79, -528)`}
        >
          <use href="#mapmarker" />
        </g>
      </Marker>
    );
  }
}

GWUMarker.propTypes = {
  marker: PropTypes.shape({
    name: PropTypes.string.isRequired,
    // [lng, lat]
    coordinates: PropTypes.arrayOf(PropTypes.number).isRequired
  }).isRequired,
  scale: PropTypes.number.isRequired
};

GWUMarker.defaultProps = {
  scale: 0.09
};

export default GWUMarker;
