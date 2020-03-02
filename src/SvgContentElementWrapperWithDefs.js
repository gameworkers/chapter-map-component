import React from "react";
import PropTypes from "prop-types";

const red = "#e5261a";

function GWUMarkerPath({ fill, stroke, strokeWidth }) {
  return (
    <path
      fill={fill}
      stroke={stroke}
      strokeWidth={strokeWidth}
      d={`
        M-78.736,56.784c-89.734,0-162.725,73.003-162.725,162.732c0,36.571,25.572,94.673,76.01,172.686
        c36.714,56.783,73.96,104.287,74.33,104.75c2.985,3.807,7.553,6.018,12.379,6.018c4.829,0,9.396-2.211,12.381-6.018
        c0.37-0.463,37.622-47.971,74.33-104.75c50.438-78.013,76.012-136.121,76.012-172.686
        C83.993,129.787,10.999,56.784-78.736,56.784z
        M-140.173,241.81c-3.481,2.971-7.626,5.503-12.439,7.6
        c-4.804,2.094-10.303,3.146-16.504,3.146c-5.293,0-10.114-0.871-14.47-2.616c-4.357-1.741-8.101-4.16-11.226-7.259
        c-3.136-3.098-5.571-6.775-7.311-11.037c-1.742-4.259-2.616-8.903-2.616-13.936v-0.197c0-4.838,0.89-9.372,2.662-13.602
        c1.772-4.222,4.248-7.937,7.405-11.131c3.16-3.193,6.909-5.712,11.234-7.55c4.317-1.839,9.061-2.758,14.228-2.758
        c3.033,0,5.791,0.208,8.28,0.631c2.481,0.417,4.787,1.016,6.916,1.788c2.136,0.776,4.128,1.745,6.002,2.907
        c1.876,1.158,3.682,2.481,5.421,3.966l-9.391,11.327c-1.291-1.098-2.583-2.067-3.868-2.906c-1.293-0.834-2.63-1.548-4.018-2.126
        c-1.395-0.583-2.891-1.033-4.507-1.356c-1.614-0.322-3.387-0.484-5.316-0.484c-2.717,0-5.247,0.565-7.603,1.697
        c-2.354,1.13-4.419,2.651-6.2,4.562c-1.771,1.906-3.158,4.139-4.158,6.692c-1.001,2.555-1.498,5.289-1.498,8.2v0.193
        c0,3.104,0.497,5.968,1.498,8.588c1,2.617,2.418,4.898,4.263,6.838c1.835,1.941,4,3.44,6.481,4.513
        c2.481,1.068,5.247,1.598,8.28,1.598c5.545,0,10.226-1.366,14.03-4.096v-9.749h-14.999v-12.872h29.423V241.81z
        M-51.398,251.879
        h-12.977l-15.29-44.337l-15.299,44.337h-12.974l-23.13-68.245h15.969l14.038,45.883l15.195-46.077h12.779l15.195,46.077
        l14.039-45.883h15.58L-51.398,251.879z
        M38.676,235.323c-1.395,3.81-3.381,6.988-5.956,9.536
        c-2.585,2.553-5.712,4.455-9.39,5.711c-3.68,1.26-7.777,1.887-12.298,1.887c-9.099,0-16.262-2.517-21.49-7.551
        c-5.224-5.033-7.839-12.55-7.839-22.554v-38.718h14.905v38.333c0,5.549,1.291,9.729,3.876,12.532
        c2.575,2.809,6.161,4.217,10.744,4.217c4.579,0,8.161-1.359,10.746-4.07c2.575-2.71,3.868-6.775,3.868-12.194v-38.817h14.912
        v38.234C40.755,227.032,40.061,231.518,38.676,235.323z
      `}
    />
  );
}

function UnionMarkerPath({ fill, stroke, strokeWidth }) {
  return (
    <path
      fill={fill}
      stroke={stroke}
      strokeWidth={strokeWidth}
      d={`
        M-78.736,56.784c-89.734,0-162.725,73.003-162.725,162.732c0,36.571,25.572,94.673,76.01,172.686
        c36.714,56.783,73.96,104.287,74.33,104.75c2.985,3.807,7.553,6.018,12.379,6.018c4.829,0,9.396-2.211,12.381-6.018
        c0.37-0.463,37.622-47.971,74.33-104.75c50.438-78.013,76.012-136.121,76.012-172.686
        C83.993,129.787,10.999,56.784-78.736,56.784z
        M-140.173,241.81c-3.481,2.971-7.626,5.503-12.439,7.6     
      `}
    />
  );
}

function SvgContentElementWrapperWithDefs({
  children,
  forceGrayscale,
  ...rest
}) {
  return (
    <React.Fragment>
      <defs>
        <filter id="grayscale">
          <feColorMatrix type="saturate" values="0" />
        </filter>
        <pattern
          id="redpattern"
          patternUnits="userSpaceOnUse"
          width="8"
          height="8"
        >
          <rect style={{ opacity: 0.7 }} width="8" height="8" fill={red} />
        </pattern>
        <pattern
          id="lessredpattern"
          patternUnits="userSpaceOnUse"
          width="8"
          height="8"
        >
          <rect style={{ opacity: 0.3 }} width="8" height="8" fill={red} />
        </pattern>
        <pattern
          id="hardlyredpattern"
          patternUnits="userSpaceOnUse"
          width="8"
          height="8"
        >
          <rect style={{ opacity: 0.04 }} width="8" height="8" fill={red} />
        </pattern>
        <g id="gwuMarker">
          <rect
            x="-218.979"
            y="163.179"
            fill="#fff"
            width="279.333"
            height="106.667"
          />
          <g>
            <GWUMarkerPath stroke="#fff" strokeWidth={20} />
            <GWUMarkerPath />
          </g>
        </g>
        <g id="unionMarker">
          <rect
            x="-218.979"
            y="163.179"
            fill="#000"
            width="279.333"
            height="106.667"
          />
          <g>
            <UnionMarkerPath fill="#000" stroke="#000" strokeWidth={20} />
            <UnionMarkerPath fill="#fff" stroke="#000" />
          </g>
        </g>
      </defs>
      <g filter={forceGrayscale ? "url(#grayscale)" : undefined}>
        {React.cloneElement(children, rest)}
      </g>
    </React.Fragment>
  );
}

SvgContentElementWrapperWithDefs.propTypes = {
  children: PropTypes.element.isRequired,
  forceGrayscale: PropTypes.bool
};

export default SvgContentElementWrapperWithDefs;
