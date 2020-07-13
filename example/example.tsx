import React, { useState, useCallback } from "react";
import { render } from "react-dom";

import ChapterMap, { ChapterMapProps } from "../src";

import DownloadButtons from "./map-download";

import "react-popper-tooltip/dist/styles.css";

import "./example.css";
import "./mapzoombuttons.css";

const MapWithDownloadButtons = (props: ChapterMapProps) => {
  const [chapterMap, setChapterMap] = useState<HTMLDivElement | null>(null);

  return (
    <div style={{ position: "relative", marginBottom: 30 }}>
      <DownloadButtons map={chapterMap!} />
      <ChapterMap {...props} ref={setChapterMap} className="chapter_map" />
    </div>
  );
};

const MapWithZoomButtons = () => {
  const [
    {
      coordinates: [x, y],
      zoom,
    },
    setPos,
  ] = useState({ coordinates: [13, 18], zoom: 0.9 });

  const handleZoomIn = useCallback(() => {
    setPos(({ coordinates, zoom: z }) => ({
      coordinates,
      zoom: Math.min(z + 0.5, 20),
    }));
  }, [setPos]);

  const handleZoomOut = useCallback(() => {
    setPos(({ coordinates, zoom: z }) => ({
      coordinates,
      zoom: Math.max(z - 0.5, 1),
    }));
  }, [setPos]);

  return (
    <div className="zoomcontainer">
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
        x={x}
        y={y}
        zoom={zoom}
        onPanZoom={setPos}
        panZoomControls
        height={450}
        markerScale={0.1}
        scale={205}
        width={780}
        tooltipClassName="zoom_tooltip"
      />
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
      <MapWithDownloadButtons
        x={15}
        y={15}
        height={551}
        markerScale={0.075}
        scale={205}
        width={980}
        zoom={1}
        forceGrayscale={grayscale}
      />
      <MapWithDownloadButtons
        x={-80}
        y={13}
        width={600}
        height={1000}
        scale={400}
        markerScale={0.08}
        geographyFilter={(geo) => geo.properties.REGION_UN === "Americas"}
        forceGrayscale={grayscale}
      />
      <MapWithDownloadButtons
        x={14}
        y={52}
        height={725}
        width={720}
        markerScale={0.09}
        scale={1125}
        zoom={0.95}
        geographyFilter={(geo) => geo.properties.REGION_UN === "Europe"}
        forceGrayscale={grayscale}
      />
      <div style={{ marginBottom: 30 }}>
        <p>Pan and zoom by dragging on the map (uncontrolled component)</p>
        <ChapterMap
          x={15}
          y={15}
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
      <div style={{ marginBottom: 30 }}>
        <p>
          Pan and zoom by dragging on the map or using buttons (controlled
          component)
        </p>
        <MapWithZoomButtons />
      </div>
    </>
  );
};

render(<DemoApp />, document.querySelector("#app"));
