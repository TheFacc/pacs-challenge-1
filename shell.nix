{ pkgs ? import <nixpkgs> {} }:

  pkgs.mkShell {
      nativeBuildInputs = with pkgs.buildPackages; [
        gcc11
        cmake
        muparser
        #muparserx
        #eigen
      ];
  }
