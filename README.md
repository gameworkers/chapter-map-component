# GWU Chapter Map Component

## See examples [here](https://gameworkers.github.io/chapter-map-component)!

Install by doing:

```console
npm install @gameworkers/chapter-map-component
# or yarn add @gameworkers/chapter-map-component
```

```js
import ChapterMapComponent from '@gameworkers/chapter-map-component';

// later...

<ChapterMapComponent
  centerLat={55}
  centerLng={15}
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
{
  mapDataUrl?: string;
  memberDataUrl?: string;
  centerLat?: number;
  centerLng?: number;
  width?: number;
  height?: number;
  scale?: number;
  geographyFilter?: (geography) => boolean;
  markerScale?: number;
  panZoomControls?: boolean;
  forceGrayscale?: boolean;
  className?: string;
  tooltipClassName?: string;
  zoom?: number;
  projection?: string;
}
```

Props have the following default values:

```js
{
  mapDataUrl = DEFAULT_MAP_DATA_URL,
  memberDataUrl = DEFAULT_MEMBER_DATA_URL,
  centerLat = 0,
  centerLng = 0,
  width = 980,
  height = 551,
  scale = 205,
  geographyFilter = (geo) => geo.properties.REGION_UN !== "Antarctica",
  markerScale = 0.09,
  panZoomControls = false;
  forceGrayscale = false,
  tooltipClassName = "gwu_chapter_tooltip",
  zoom = 1,
  projection = "geoNaturalEarth1"
}
```
