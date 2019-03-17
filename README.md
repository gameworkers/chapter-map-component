# GWU Chapter Map Component

## See examples [here](https://gameworkers.github.io/chapter-map-component)!

Install by doing:

```console
npm install --save @gameworkers/chapter-map-component
```

```js
import ChapterMapComponent from '@gameworkers/chapter-map-component';

window.WORLD_110M_JSON_PATH = location.pathname + 'dist/world-100m.json';
window.WORLD_50M_JSON_PATH = location.pathname + 'dist/world-50m.json';

// later...

<ChapterMapComponent
  centerLat={55}
  centerLng={15}
  className="chapter_map"
  forceGrayscale={false}
  height={825}
  isGeographyIncluded={function(geography) {
    return geography.properties.REGION_UN === 'Europe';
  }}
  markerScale={0.1}
  scale={1125}
  width={720}
/>
```

`ChapterMapComponent` expects [`world-110m.json`](/src/world-110m.json) and [`world-50m.json`](/src/world-50m.json) to be served from your application's web server, and paths declared in the global scope on the client side via `WORLD_110M_JSON_PATH` and `WORLD_50M_JSON_PATH`.

All props are optional and have the following default values if nothing is specified:

```js
ChapterMapComponent.defaultProps = {
  centerLat: 0,
  centerLng: 0,
  width: 980,
  height: 551,
  scale: 205,
  isGeographyIncluded: () => true,
  markerScale: 0.1,
  forceGrayscale: false
};
```
