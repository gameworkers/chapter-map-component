# GWU Chapter Map Component

**Check out examples [here](https://gameworkers.github.io/chapter-map-component)!**

Install by doing:

```console
npm install @gameworkers/chapter-map-component
# or yarn add @gameworkers/chapter-map-component
```

```js
import ChapterMapComponent from '@gameworkers/chapter-map-component';

// import the tooltip styles via a method of your choosing.
import 'react-popper-tooltip/dist/styles.css';

// later...

<ChapterMapComponent
  x={15}
  y={55}
  className="chapter_map"
  forceGrayscale={false}
  height={825}
  geographyFilter={
    (geo) => geo.properties.REGION_UN === 'Europe'
  }
  markerScale={0.1}
  scale={1125}
  width={720}
/>
```

The component accepts the following props:

```ts
interface ChapterMapProps {
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
```

Props have the following default values:

```js
{
  x = 0,
  y = 0,
  zoom = 1,
  panZoomControls = false,
  mapDataUrl = DEFAULT_MAP_DATA_URL,
  memberDataUrl = DEFAULT_MEMBER_DATA_URL,
  width = 980,
  height = 551,
  scale = 205,
  geographyFilter = (geo) => geo.properties.REGION_UN !== "Antarctica",,
  markerScale = 0.09,
  forceGrayscale = false,
  tooltipClassName = "gwu_chapter_tooltip",
  projection = "geoNaturalEarth1"
}
```
