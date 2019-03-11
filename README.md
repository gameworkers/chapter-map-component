# GWU Chapter Map Component

## See examples [here](https://gameworkers.github.io/chapter-map-component)!

Install by doing:

```console
npm install --save @gameworkers/chapter-map-component
```

```js
import ChapterMapComponent from '@gameworkers/chapter-map-component';

// later...

<ChapterMapComponent
  centerLat={55}
  centerLng={15}
  className="chapter_map"
  height={825}
  isGeographyIncluded={function(geography) {
    return geography.properties.REGION_UN === 'Europe';
  }}
  markerScale={0.1}
  scale={1125}
  width={720}
/>
```
