import React, { memo, forwardRef, ForwardRefRenderFunction } from "react";
import { Marker } from "@sux/react-simple-maps";

interface ExtraEvents {
  onTouchEnd?: (event: React.SyntheticEvent) => void;
  onClick?: (event: React.SyntheticEvent) => void;
  onContextMenu?: (event: React.SyntheticEvent) => void;
  onMouseEnter?: (event: React.SyntheticEvent) => void;
  onMouseLeave?: (event: React.SyntheticEvent) => void;
  onMouseMove?: (event: React.SyntheticEvent) => void;
  onFocus?: (event: React.SyntheticEvent) => void;
  onBlur?: (event: React.SyntheticEvent) => void;
}

type GWUMarkerProps = {
  marker: {
    name: string;
    coordinates: [number, number];
    data: {
      isUnion?: boolean;
    };
  };
  scale?: number;
} & ExtraEvents;

const GWUMarker: ForwardRefRenderFunction<SVGGElement, GWUMarkerProps> = (
  { marker, scale = 0.09, ...markerProps },
  ref
) => (
  <Marker key={marker.name} coordinates={marker.coordinates} {...markerProps}>
    <g
      ref={ref}
      style={{ cursor: "pointer" }}
      transform={`scale(${scale}), translate(79, -528)`}
    >
      <use href={marker.data.isUnion ? "#unionMarker" : "#gwuMarker"} />
    </g>
  </Marker>
);

export default memo(forwardRef(GWUMarker));
