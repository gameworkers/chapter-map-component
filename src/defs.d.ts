declare module "@sux/react-simple-maps" {
  export * from "react-simple-maps";
  import { ZoomableGroupProps } from "react-simple-maps";

  // these are on the normal ZoomableGroupProps too, but we can override the
  // incorrect typings more easily like this
  export type CustomZoomableGroupProps = Omit<
    ZoomableGroupProps,
    "onMoveStart" | "onMoveEnd" | "onMove"
  > & {
    onMoveStart?: (
      pos: { coordinates: [number, number]; zoom: number },
      d3Event: Event
    ) => void;
    onMove?: (
      pos: { x: number; y: number; k: number; dragging: React.DragEvent },
      d3Event: Event
    ) => void;
    onMoveEnd?: (
      pos: { coordinates: [number, number]; zoom: number },
      d3Event: Event
    ) => void;
    filterZoomEvent?: (d3Event: Event) => boolean;
  };

  export const CustomZoomableGroup: React.FunctionComponent<CustomZoomableGroupProps>;
}
