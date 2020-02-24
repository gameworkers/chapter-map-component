var CompressionPlugin = require('compression-webpack-plugin');
var path = require('path');

var OUTPUT_DIR = './dist';

function webpackConfig({
  fileName,
  entryPath,
  minimize,
  externals = {},
  libraryName = undefined
}) {
  const publicPath = '/dist';
  const config = {
    mode: process.env.WEBPACK_SERVE ? 'development' : 'production',
    entry: {
      [minimize ? `${fileName}.min` : fileName]: entryPath
    },
    resolve: {
      extensions: ['.js', '.jsx']
    },
    output: {
      path: path.join(process.cwd(), OUTPUT_DIR),
      publicPath,
      libraryTarget: 'umd',
      libraryExport: 'default',
      library: libraryName,
      filename: '[name].js',
      globalObject: '(typeof self !== "undefined" ? self : this)'
    },
    serve: {
      devMiddleware: {
        publicPath
      }
    },
    module: {
      rules: [
        {
          test: /\.jsx?$/,
          exclude: /node_modules/,
          use: { loader: 'babel-loader' }
        }
      ]
    },
    externals: {
      'prop-types': {
        root: 'PropTypes',
        commonjs: 'prop-types',
        commonjs2: 'prop-types',
        amd: 'prop-types'
      },
      react: {
        root: 'React',
        commonjs: 'react',
        commonjs2: 'react',
        amd: 'react'
      },
      ...externals
    },
    devtool: 'source-map',
    plugins: [
      new CompressionPlugin(),
    ],
    optimization: {
      noEmitOnErrors: true,
      minimize
    }
  };
  return config;
}

let configs = [
  {
    fileName: 'chapter-map-component',
    entryPath: './src/index.js',
    minimize: false,
    libraryName: 'ChapterMapComponent',
    externals: process.env.WEBPACK_SERVE
      ? {}
      : {
        '@gameworkers/react-simple-maps': {
          root: 'ReactSimpleMaps',
          commonjs: '@gameworkers/react-simple-maps',
          commonjs2: '@gameworkers/react-simple-maps',
          amd: '@gameworkers/react-simple-maps'
        },
        'redux': {
          root: 'Redux',
          commonjs: 'redux',
          commonjs2: 'redux',
          amd: 'redux'
        },
        'next-redux-wrapper': {
          root: 'NextReduxWrapper',
          commonjs: 'next-redux-wrapper',
          commonjs2: 'next-redux-wrapper',
          amd: 'next-redux-wrapper'
        },
        'redux-tooltip': {
          root: 'ReduxTooltip',
          commonjs: 'redux-tooltip',
          commonjs2: 'redux-tooltip',
          amd: 'redux-tooltip'
        },
        'react-redux': {
          root: 'ReactRedux',
          commonjs: 'react-redux',
          commonjs2: 'react-redux',
          amd: 'react-redux'
        }
      }
  },
  {
    fileName: 'demo-externals',
    entryPath: './demo-externals.js',
    minimize: false
  }
];
if (process.env.BUILD_MODE !== 'unminimized') {
  configs = configs.concat(configs.map(c => ({ ...c, minimize: true })));
}

module.exports = configs.map(webpackConfig);
