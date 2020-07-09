import React, {
  useState,
  useEffect,
  forwardRef,
  memo,
  ForwardRefRenderFunction,
} from "react";

import TooltipTrigger, { TooltipArg, ChildrenArg } from "react-popper-tooltip";

import {
  ComposableMap,
  ZoomableGroup,
  Geographies,
  Geography,
} from "react-simple-maps";
// } from "@gameworkers/react-simple-maps";

import useFetch from "./use-fetch";
import SvgContentElementWrapperWithDefs from "./SvgContentElementWrapperWithDefs";
import GWUMarker from "./GWUMarker";

const DEFAULT_MAP_DATA_URL =
  "https://gameworkers.github.io/data/third_party/world-50m.json";
const DEFAULT_MEMBER_DATA_URL =
  "https://gameworkers.github.io/data/members.json";

interface Member {
  isUnion?: boolean;
  isChapter?: boolean;
  lat: number;
  lng: number;
  location: string;
  country?: string;
  chapterInfo: {
    description?: string;
    applicationLink?: string;
    twitter?: string;
    email?: string;
    website?: string;
  };
}
interface Marker {
  name: string;
  coordinates: [number, number];
  data: Member;
}

interface Geo {
  properties: { [key: string]: string };
}

const countryNameKeys = [
  "ABBREV",
  "CONTINENT",
  "FORMAL_EN",
  "ISO_A2",
  "ISO_A3",
  "NAME",
  "NAME_LONG",
] as const;

const geographyMatchesCountryString = (
  geography: Geo,
  countryString?: string
) => countryNameKeys.some((key) => geography.properties[key] === countryString);

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
  /**
   * URL where map data JSON will be fetched from. The URL should end with a
   * slash and will be concatenated with the required filenames (eg.
   * "world-50m.json", "world-110m.json", "members.json") to obtain the final
   * fetch URL.
   */
  mapDataUrl?: string;
  memberDataUrl?: string;
  centerLat?: number;
  centerLng?: number;
  width?: number;
  height?: number;
  scale?: number;
  geographyFilter?: (geography: Geo) => boolean;
  markerScale?: number;
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
          <ZoomableGroup center={[centerLng, centerLat]} zoom={zoom}>
            <Geographies geography={mapData}>
              {({ geographies }: { geographies: Geo[] }) =>
                geographies.filter(geographyFilter).map((geography, i) => {
                  let hasMatchingPoint = false;
                  if (memberData) {
                    hasMatchingPoint = memberData.some((member) => {
                      return geographyMatchesCountryString(
                        geography,
                        member.country
                      );
                    });
                  }

                  const style = {
                    fill: hasMatchingPoint
                      ? "url(#redpattern)"
                      : "url(#hardlyredpattern)",
                    stroke: "#222",
                    strokeWidth: 0.5 / zoom,
                    outline: "none",
                  };

                  return (
                    <Geography
                      key={i}
                      geography={geography}
                      style={{
                        default: style,
                        hover: style,
                        pressed: style,
                      }}
                    />
                  );
                })
              }
            </Geographies>

            {markers?.map((marker) => (
              <TooltipTrigger
                key={marker.name}
                placement="top"
                trigger="click"
                tooltip={(args: TooltipArg) => (
                  <Tooltip {...args}>
                    {getTooltipContent(marker, tooltipClassName)}
                  </Tooltip>
                )}
              >
                {({ triggerRef, getTriggerProps }: ChildrenArg) => (
                  <GWUMarker
                    ref={triggerRef as any}
                    marker={marker}
                    scale={markerScale}
                    {...getTriggerProps()}
                  />
                )}
              </TooltipTrigger>
            ))}
          </ZoomableGroup>
        </SvgContentElementWrapperWithDefs>
      </ComposableMap>
    </div>
  );
};

const Tooltip = ({
  arrowRef,
  tooltipRef,
  getArrowProps,
  getTooltipProps,
  placement,
  children,
}: TooltipArg & { children: React.ReactNode }) => (
  <div
    {...getTooltipProps({
      ref: tooltipRef,
      className: "tooltip-container",
    })}
  >
    <div
      {...getArrowProps({
        ref: arrowRef,
        className: "tooltip-arrow",
        "data-placement": placement,
      })}
    />
    {children}
  </div>
);

const getTooltipContent = (marker: Marker, className?: string) => {
  const {
    data: { location, chapterInfo },
  } = marker;
  if (chapterInfo) {
    const {
      description,
      applicationLink,
      twitter,
      email,
      website,
    } = chapterInfo;

    return (
      <div className={className}>
        <h3>{location}</h3>
        {description && <p>{description}</p>}
        {(twitter || email || website) && (
          <p>
            {twitter && (
              <>
                <strong>Twitter: </strong>
                <a href={`https://twitter.com/${twitter}`}>@{twitter}</a>
                <br />
              </>
            )}
            {email && (
              <>
                <strong>Email: </strong>
                <a href={`mailto:${email}`}>{email}</a>
                <br />
              </>
            )}
            {website && (
              <>
                <strong>Website: </strong>
                <a
                  href={`${
                    website.startsWith("http") ? "" : "http://"
                  }${website}`}
                >
                  {website}
                </a>
                <br />
              </>
            )}
          </p>
        )}
        {applicationLink && (
          <p>
            <a href={applicationLink}>Apply here!</a>
          </p>
        )}
      </div>
    );
  }
};

export default memo(forwardRef(ChapterMap));
