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

// bad typings
declare module "react-simple-maps" {
  // eslint-disable-next-line no-shadow
  export const CustomZoomableGroup: React.FunctionComponent<ZoomableGroupProps>;
  interface ZoomableGroupProps {
    filterZoomEvent: (d3Event: Event) => boolean;
  }
}

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
  mapDataUrl?: string;
  memberDataUrl?: string;
  centerLat?: number;
  centerLng?: number;
  width?: number;
  height?: number;
  scale?: number;
  geographyFilter?: (geography: Geo) => boolean;
  markerScale?: number;
  panZoomControls?: boolean;
  forceGrayscale?: boolean;
  className?: string;
  tooltipClassName?: string;
  zoom?: number;
  projection?: string;
}

const ChapterMap: ForwardRefRenderFunction<HTMLDivElement, ChapterMapProps> = (
  {
    mapDataUrl = DEFAULT_MAP_DATA_URL,
    memberDataUrl = DEFAULT_MEMBER_DATA_URL,
    centerLat = 0,
    centerLng = 0,
    width = 980,
    height = 551,
    scale = 205,
    geographyFilter = defaultGeographyFilter,
    markerScale = 0.09,
    panZoomControls = false,
    forceGrayscale = false,
    tooltipClassName = "gwu_chapter_tooltip",
    zoom = 1,
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
            center={[centerLng, centerLat]}
            zoom={zoom}
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
