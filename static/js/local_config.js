(function() {
  citynavi.update_configs({
    defaults: {
      new_feature_api_url: "http://example.com/"
    },
    helsinki: {
      cities: citynavi.configs.helsinki.cities.concat(["Inarinjärvi"])
    },
    llanfairpwll: {
      country: "gb"
    },
    overrides: {
      osm_notes_url: "http://api06.dev.openstreetmap.org/api/0.6/notes.json"
    }
  });

  citynavi.set_config("helsinki");

}).call(this);

//# sourceMappingURL=local_config.js.map
