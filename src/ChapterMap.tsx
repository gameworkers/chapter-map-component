import React, {
  useState,
  useEffect,
  forwardRef,
  memo,
  ForwardRefRenderFunction,
} from "react";

import TooltipTrigger, { TooltipArg, ChildrenArg } from "react-popper-tooltip";
import useSWR from "swr";

import {
  ComposableMap,
  ZoomableGroup,
  Geographies,
  Geography,
} from "react-simple-maps";
// } from "@gameworkers/react-simple-maps";

import SvgContentElementWrapperWithDefs from "./SvgContentElementWrapperWithDefs";
import GWUMarker from "./GWUMarker";

const DEFAULT_MAP_DATA_URL = "https://gameworkers.github.io/data/";

const WORLD_DATA = "third_party/world-50m.json";
const WORLD_DATA_SMALL = "third_party/world-110m.json";
const MEMBER_DATA = "members.json";

const wrapperStyles = {
  width: "100%",
  maxWidth: 980,
  margin: "0 auto",
  fontFamily: "Roboto, sans-serif",
};

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
    .map((member) => ({
      name: member.location,
      coordinates: [member.lng, member.lat],
      data: member,
    }));

// const getContentForTooltip = (marker: Marker) => {
//   const {
//     data: { location, chapterInfo },
//   } = marker;
//   let content: string = location;
//   if (chapterInfo) {
//     const {
//       description,
//       applicationLink,
//       twitter,
//       email,
//       website,
//     } = chapterInfo;
//     content = `<h3>${location}</h3>`;
//     if (description) {
//       content += `<p>${description}</p>`;
//     }
//     if (twitter || email || website) {
//       content += "<p>";
//     }
//     if (twitter) {
//       content += `
//       <strong>Twitter:</strong>
//       <a href="https://twitter.com/${twitter}">
//         @${twitter}
//       </a>
//       <br />
//     `;
//     }
//     if (email) {
//       content += `
//       <strong>Email:</strong>
//       <a href="mailto:${email}">
//         ${email}
//       </a>
//       <br />
//     `;
//     }
//     if (website) {
//       content += `
//       <strong>Website:</strong>
//       <a href="${website.indexOf("http") === 0 ? "" : "http://"}${website}">
//         ${website}
//       </a>
//       <br />
//     `;
//     }
//     if (twitter || email || website) {
//       content += "</p>";
//     }
//     if (applicationLink) {
//       content += `
//       <p>
//         <a href="${applicationLink}">
//           Apply here!
//         </a>
//       </p>
//     `;
//     }
//     content = `<div class="${tooltipClassName}">${content}</div>`;
//   }
//   dispatchTooltip(evt, content);
// };

interface ChapterMapProps {
  /**
   * URL where map data JSON will be fetched from. The URL should end with a
   * slash and will be concatenated with the required filenames (eg.
   * "world-50m.json", "world-110m.json", "members.json") to obtain the final
   * fetch URL.
   */
  mapDataUrl?: string;
  centerLat?: number;
  centerLng?: number;
  width?: number;
  height?: number;
  scale?: number;
  isGeographyIncluded?: (geography: Geo) => boolean;
  markerScale?: number;
  forceGrayscale?: boolean;
  className?: string;
  tooltipClassName?: string;
  zoom?: number;
  enablePanning?: boolean;
  projection?: string;
}

const ChapterMap: ForwardRefRenderFunction<HTMLDivElement, ChapterMapProps> = (
  {
    mapDataUrl = DEFAULT_MAP_DATA_URL,
    centerLat = 0,
    centerLng = 0,
    width = 980,
    height = 551,
    scale = 205,
    isGeographyIncluded = (geography) =>
      geography.properties.REGION_UN !== "Antarctica",
    markerScale = 0.09,
    forceGrayscale = false,
    tooltipClassName = "gwu_chapter_tooltip",
    zoom = 1,
    enablePanning = false,
    projection = "geoNaturalEarth1",
    className,
  },
  ref
) => {
  const useFetch = <T,>(key: string) =>
    useSWR<T>(
      mapDataUrl + key,
      (url: string) => fetch(url).then((res) => res.json()),
      {
        revalidateOnFocus: false,
        revalidateOnReconnect: false,
        refreshWhenOffline: false,
        refreshWhenHidden: false,
        refreshInterval: 0,
      }
    );

  const results = [
    useFetch<{ [k: string]: any }>(WORLD_DATA_SMALL),
    useFetch<{ [k: string]: any }>(WORLD_DATA),
    useFetch<Member[]>(MEMBER_DATA),
  ] as const;

  const [
    { data: smallMapData },
    { data: mapData },
    { data: memberData },
  ] = results;

  const errors = results
    .map((r) => r.error)
    .filter((e) => e)
    .join(", ");

  const bestAvailableMapData = mapData || smallMapData;

  const [markers, setMarkers] = useState<Marker[] | null>(null);

  // const focusedMarker = useRef<HTMLElement | null>(null);

  useEffect(() => {
    if (!memberData) {
      setMarkers(null);
      return;
    }
    setMarkers(getMarkerData(memberData));
  }, [memberData]);

  // useEffect(() => {
  //   window.addEventListener("click", handleOutsideMarkerClick);
  //   return () => window.removeEventListener("click", handleOutsideMarkerClick);
  // });

  // const dispatchTooltip = useCallback((evt: React.MouseEvent, content) => {
  //   const x = evt.clientX;
  //   const y = evt.clientY + window.pageYOffset;
  //   dispatch(
  //     actions.show({
  //       origin: { x, y },
  //       content,
  //     })
  //   );
  // }, []);

  // const hideTooltip = useCallback(() => {
  //   dispatch(actions.hide());
  //   Promise.resolve().then(() => {
  //     focusedMarker.current = null;
  //   });
  // }, []);

  // const handleCountryMove = useCallback(
  //   (geography, evt) => {
  //     if (!focusedMarker.current) {
  //       dispatchTooltip(evt, geography.properties.NAME_LONG);
  //     }
  //   },
  //   [dispatchTooltip]
  // );

  // const handleCountryLeave = useCallback(() => {
  //   if (focusedMarker.current) hideTooltip();
  // }, [hideTooltip]);

  // const handleMarkerMove = useCallback(
  //   (marker, evt) => {
  //     if (!focusedMarker.current) dispatchTooltip(evt, marker.name);
  //   },
  //   [dispatchTooltip]
  // );

  // const handleMarkerLeave = useCallback(() => {
  //   if (focusedMarker.current) hideTooltip();
  // }, [hideTooltip]);

  // const handleMarkerClick = useCallback(
  //   (marker: Marker, evt: React.MouseEvent) => {
  //     evt.stopPropagation();
  //     evt.nativeEvent.stopImmediatePropagation();

  //     focusedMarker.current = focusedMarker.current === marker ? null : marker;

  //     if (focusedMarker.current) {
  //     } else {
  //       hideTooltip();
  //     }
  //   },
  //   [dispatchTooltip, hideTooltip, tooltipClassName]
  // );

  // const handleOutsideMarkerClick = hideTooltip;

  if (errors) {
    console.error(errors);
    return <div>{errors}</div>;
  }

  if (!bestAvailableMapData) {
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
            <Geographies geography={bestAvailableMapData}>
              {({ geographies }: { geographies: Geo[] }) =>
                geographies.filter(isGeographyIncluded).map((geography, i) => {
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
                // usePortal={false}
                tooltip={Tooltip}
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
}: TooltipArg) => (
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
    Hello, World!
  </div>
);

export default memo(forwardRef(ChapterMap));
