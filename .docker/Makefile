# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>

SHELL = /bin/bash
.ONESHELL:

.PHONY: default build push login

default:
	@echo Valid Makefile targets are build and push

build:
	@docker image inspect platipodium/tame > /dev/null 2>&1 || docker build -f ./Dockerfile-fabm -t platipodium/tame ..

login:
	docker login -u ${CI_REGISTRY_USER} -p ${CI_WRITE_REGISTRY_TOKEN} ${CI_REGISTRY}

push: build login
	docker image tag $(shell docker images -q platipodium/tame) ${CI_REGISTRY}/kse/tame:tame && \
	docker image push ${CI_REGISTRY}/kse/tame:tame
