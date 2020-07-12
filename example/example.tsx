import React, { useState, useRef, useCallback } from "react";
import { render } from "react-dom";

import ChapterMap, { ChapterMapProps } from "../src";

import DownloadButtons from "./map-download";

import "./example.css";
import "./mapzoombuttons.css";

import "../tooltip-styles.css";

const MapWithDownloadButtons = (props: ChapterMapProps) => {
  const chapterMap = useRef<HTMLDivElement | null>(null);

  return (
    <div style={{ position: "relative", marginBottom: 15 }}>
      <DownloadButtons map={chapterMap.current!} />
      <ChapterMap {...props} ref={chapterMap} className="chapter_map" />
    </div>
  );
};

const MapWithZoomButtons = () => {
  const [zoom, setZoom] = useState(0.9);

  const handleZoomIn = useCallback(() => {
    setZoom((z) => Math.min(z + 0.5, 20));
  }, [setZoom]);

  const handleZoomOut = useCallback(() => {
    setZoom((z) => Math.max(z - 0.5, 1));
  }, [setZoom]);

  return (
    <div className="zoomcontainer" style={{ marginBottom: 30 }}>
      <div className="zoom_buttonscontainer">
        <div className="zoom_buttons">
          <button
            onClick={handleZoomIn}
            disabled={zoom >= 20}
            style={{
              marginBottom: 10,
              opacity: zoom >= 20 ? 0.3 : 1,
            }}
            className="zoom_button"
          >
            ➕
          </button>
          <button
            onClick={handleZoomOut}
            disabled={zoom <= 1}
            style={{
              opacity: zoom <= 1 ? 0.3 : 1,
            }}
            className="zoom_button"
          >
            ➖
          </button>
        </div>
      </div>
      <ChapterMap
        className="chapter_map"
        centerLat={18}
        centerLng={13}
        height={450}
        markerScale={0.1}
        scale={205}
        width={780}
        panZoomControls
        zoom={zoom}
        tooltipClassName="zoom_tooltip"
      />
    </div>
  );
};

const DemoApp = () => {
  const [grayscale, setGrayscale] = useState(false);

  return (
    <>
      <MapWithZoomButtons />

      <button
        style={{ marginBottom: 20 }}
        onClick={() => setGrayscale(!grayscale)}
      >
        Toggle Grayscale
      </button>
      <MapWithDownloadButtons
        centerLat={15}
        centerLng={15}
        height={551}
        markerScale={0.075}
        scale={205}
        width={980}
        zoom={1}
        forceGrayscale={grayscale}
      />
      <MapWithDownloadButtons
        centerLat={13}
        centerLng={-80}
        width={600}
        height={1000}
        scale={400}
        markerScale={0.08}
        geographyFilter={(geo) => geo.properties.REGION_UN === "Americas"}
        forceGrayscale={grayscale}
      />
      <MapWithDownloadButtons
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
      <div style={{ marginBottom: 15 }}>
        <p>Pannable and zoomable version:</p>
        <ChapterMap
          centerLng={15}
          centerLat={15}
          zoom={1}
          panZoomControls
          height={551}
          markerScale={0.075}
          scale={205}
          width={980}
          className="chapter_map"
          forceGrayscale={grayscale}
        />
      </div>
    </>
  );
};

render(<DemoApp />, document.querySelector("#app"));
