import React, { useState, useRef } from "react";
import { render } from "react-dom";

import ChapterMap, { ChapterMapProps } from "../src";

import DownloadButtons from "./map-download";

import "./example.css";
import "../tooltip-styles.css";

const MapWithButtons = (props: ChapterMapProps) => {
  const chapterMap = useRef<HTMLDivElement | null>(null);

  return (
    <div style={{ position: "relative", marginBottom: 15 }}>
      <DownloadButtons map={chapterMap.current!} />
      <ChapterMap {...props} ref={chapterMap} className="chapter_map" />
    </div>
  );
};

const DemoApp = () => {
  const [grayscale, setGrayscale] = useState(false);

  return (
    <>
      <button
        style={{ marginBottom: 20 }}
        onClick={() => setGrayscale(!grayscale)}
      >
        Toggle Grayscale
      </button>
      <MapWithButtons
        centerLat={15}
        centerLng={15}
        height={551}
        markerScale={0.075}
        scale={205}
        width={980}
        zoom={1}
        forceGrayscale={grayscale}
      />
      <MapWithButtons
        centerLat={13}
        centerLng={-80}
        width={600}
        height={1000}
        scale={400}
        markerScale={0.08}
        geographyFilter={(geo) => geo.properties.REGION_UN === "Americas"}
        forceGrayscale={grayscale}
      />
      <MapWithButtons
        centerLat={52}
        centerLng={14}
        height={725}
        width={720}
        markerScale={0.09}
        scale={1125}
        zoom={0.95}
        geographyFilter={(geo) => geo.properties.REGION_UN === "Europe"}
        forceGrayscale={grayscale}
      />
    </>
  );
};

render(<DemoApp />, document.querySelector("#app"));
