import React, {
  useState,
  useEffect,
  forwardRef,
  memo,
  ForwardRefRenderFunction,
} from "react";

import {
  ComposableMap,
  CustomZoomableGroup,
  Geographies,
} from "@sux/react-simple-maps";

import SvgContentElementWrapperWithDefs from "./SvgContentElementWrapperWithDefs";
import TintedGeographies from "./TintedGeographies";
import MarkerWithTooltip from "./MarkerWithTooltip";
import useFetch from "./use-fetch";

import type { Geo, Member, Marker } from "./types";

const DEFAULT_MAP_DATA_URL =
  "https://gameworkers.github.io/data/third_party/world-50m.json";
const DEFAULT_MEMBER_DATA_URL =
  "https://gameworkers.github.io/data/members.json";

const getMarkerData = (memberData: Member[]): Marker[] =>
  memberData
    .filter((member) => member.isChapter || member.isUnion)
    .sort((a, b) => {
      // we want to make sure markers lower on the map are painted in front
      if (a.lat < b.lat) return 1;
      if (a.lat > b.lat) return -1;
      if (a.lng < b.lng) return 1;
      if (a.lng > b.lng) return -1;
      return 0;
    })
    .map((member) => ({
      name: member.location,
      coordinates: [member.lng, member.lat],
      data: member,
    }));

export const defaultGeographyFilter = (geography: Geo): boolean =>
  geography.properties.REGION_UN !== "Antarctica";

export interface ChapterMapProps {
  /** The longitudinal offset of the map. */
  x?: number;
  /** The latitudinal offset of the map. */
  y?: number;
  zoom?: number;
  /**
   * Add this handler to make this a controlled component. You'll be responsible
   * for passing the x, y, and zoom back into the component;
   */
  onPanZoom?: (pos: { coordinates: [number, number]; zoom: number }) => void;
  /** Is panning/zooming permitted? */
  panZoomControls?: boolean;
  mapDataUrl?: string;
  memberDataUrl?: string;
  width?: number;
  height?: number;
  scale?: number;
  geographyFilter?: (geography: Geo) => boolean;
  markerScale?: number;
  forceGrayscale?: boolean;
  className?: string;
  tooltipClassName?: string;
  projection?: string;
}

const ChapterMap: ForwardRefRenderFunction<HTMLDivElement, ChapterMapProps> = (
  {
    x = 0,
    y = 0,
    zoom = 1,
    onPanZoom,
    panZoomControls = false,
    mapDataUrl = DEFAULT_MAP_DATA_URL,
    memberDataUrl = DEFAULT_MEMBER_DATA_URL,
    width = 980,
    height = 551,
    scale = 205,
    geographyFilter = defaultGeographyFilter,
    markerScale = 0.09,
    forceGrayscale = false,
    tooltipClassName = "gwu_chapter_tooltip",
    projection = "geoNaturalEarth1",
    className,
  },
  ref
) => {
  const { data: mapData, error: mapError } = useFetch<{ [k: string]: any }>(
    mapDataUrl
  );
  const { data: memberData, error: memberError } = useFetch<Member[]>(
    memberDataUrl
  );

  const [markers, setMarkers] = useState<Marker[] | null>(null);

  useEffect(() => {
    setMarkers(memberData ? getMarkerData(memberData) : null);
  }, [memberData]);

  if (mapError || memberError) {
    const errors = [mapError, memberError].join("\n");
    console.error("Error(s) fetching data:", errors);
    return <pre>Error: {errors}</pre>;
  }

  if (!mapData) {
    return <div>Loading...</div>;
  }

  return (
    <div ref={ref} className={className} style={{ width, height }}>
      <ComposableMap
        projection={projection}
        projectionConfig={{ scale }}
        width={width}
        height={height}
        style={{ width: "100%", height: "auto", backgroundColor: "#fff" }}
      >
        <SvgContentElementWrapperWithDefs forceGrayscale={forceGrayscale}>
          <CustomZoomableGroup
            filterZoomEvent={() => panZoomControls}
            center={[x, y]}
            zoom={zoom}
            onMoveEnd={onPanZoom}
          >
            {({ k }: { x: number; y: number; k: number }) => (
              <>
                <Geographies geography={mapData}>
                  {({ geographies }: { geographies: Geo[] }) => (
                    <TintedGeographies
                      geo={geographies}
                      zoom={k}
                      filter={geographyFilter}
                      memberData={memberData!}
                    />
                  )}
                </Geographies>
                {markers?.map((marker, i) => (
                  <MarkerWithTooltip
                    key={i}
                    marker={marker}
                    className={tooltipClassName}
                    scale={markerScale / k}
                  />
                ))}
              </>
            )}
          </CustomZoomableGroup>
        </SvgContentElementWrapperWithDefs>
      </ComposableMap>
    </div>
  );
};

export default memo(forwardRef(ChapterMap));
