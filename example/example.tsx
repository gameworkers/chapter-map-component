import React, { useState, useRef, useEffect, useCallback } from "react";
import { render } from "react-dom";

import ChapterMap from "../src";

import "../src/tooltip-styles.css";

const BIG_PNG_SIZE = 6000;

// these seem to be hardcoded to a url anyway?
// (window as any).WORLD_110M_JSON_PATH = "world-110m.json";
// (window as any).WORLD_50M_JSON_PATH = "world-50m.json";

function getSvgUrl(map: HTMLElement) {
  let content = map.querySelector("svg")!.outerHTML;
  // make sure it's a valid SVG file format!
  if (!content.includes('xmlns="http://www.w3.org/2000/svg"')) {
    content = content.replace(
      "<svg",
      '<svg xmlns="http://www.w3.org/2000/svg" version="1.1"'
    );
  }
  return URL.createObjectURL(new Blob([content], { type: "image/svg+xml" }));
}

const DownloadButtons = ({
  grayscale,
  setGrayscale,
  map,
}: {
  grayscale: boolean;
  setGrayscale: (g: boolean) => void;
  map: HTMLDivElement;
}) => {
  const anchor = useRef<HTMLAnchorElement | null>(null);
  const canvas = useRef<HTMLCanvasElement | null>(null);

  const requestPng = useCallback(
    (big?: boolean) => {
      const svgUrl = getSvgUrl(map);
      const img = new Image();

      const a = anchor.current!;
      const c = canvas.current!;
      const ctx = c.getContext("2d")!;

      img.onload = function () {
        if (big) {
          c.width = BIG_PNG_SIZE;
          c.height = (BIG_PNG_SIZE * img.height) / img.width;
        } else {
          c.width = img.width;
          c.height = img.height;
        }
        ctx.drawImage(
          img,
          0,
          0,
          img.width,
          img.height,
          0,
          0,
          c.width,
          c.height
        );
        c.toBlob((blob) => {
          const imageUrl = window.URL.createObjectURL(blob);
          a.href = imageUrl;
          a.download = `Chapter Map${big ? " (big)" : ""}.png`;
          a.click();
          URL.revokeObjectURL(svgUrl);
          URL.revokeObjectURL(imageUrl);
        });
      };
      img.src = svgUrl;
    },
    [map]
  );

  const requestSvg = useCallback(() => {
    const svgUrl = getSvgUrl(map);
    const a = anchor.current!;
    a.href = svgUrl;
    a.download = "Chapter Map.svg";
    a.click();
    URL.revokeObjectURL(svgUrl);
  }, [map]);

  useEffect(() => {
    anchor.current = document.createElement("a");
    anchor.current.style.display = "none";
    document.body.appendChild(anchor.current);

    canvas.current = document.createElement("canvas");

    return () => {
      canvas.current!.remove();
      anchor.current = null;
      canvas.current!.remove();
      canvas.current = null;
    };
  }, []);

  return (
    <div style={{ position: "absolute", top: 5, left: 5 }}>
      <button onClick={() => setGrayscale(!grayscale)}>Toggle Grayscale</button>
      <button onClick={requestSvg}>Save SVG</button>
      <button onClick={() => requestPng()}>Save PNG</button>
      <button onClick={() => requestPng(true)}>Save big PNG (for print)</button>
    </div>
  );
};

const DemoApp = () => {
  const [grayscale, setGrayscale] = useState(false);

  const chapterMap = useRef<HTMLDivElement | null>(null);

  return (
    <>
      <DownloadButtons
        grayscale={grayscale}
        setGrayscale={setGrayscale}
        map={chapterMap.current!}
      />
      <ChapterMap
        ref={chapterMap}
        centerLat={20}
        centerLng={-90}
        width={600}
        height={1000}
        scale={400}
        markerScale={0.08}
        className="chapter_map"
        forceGrayscale={grayscale}
      />
    </>
  );
};

render(<DemoApp />, document.querySelector("#app"));
