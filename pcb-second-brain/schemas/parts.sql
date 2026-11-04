-- First-cut schema for the parts store (PostgreSQL 16+, pgvector loaded separately).
-- Alembic migrations in the implementation repo are generated from this starting point.

create extension if not exists vector;
create extension if not exists pg_trgm;

create table manufacturer (
  id            bigserial primary key,
  name          text not null unique,
  aliases       text[] not null default '{}',
  website       text
);

create table package (
  id            bigserial primary key,
  name          text not null unique,          -- e.g. 'QFN-48 7x7 P0.5 EP'
  jedec_ref     text,
  ipc7351_name  text,
  pitch_mm      numeric,
  body_l_mm     numeric,
  body_w_mm     numeric,
  thermal_pad   boolean not null default false
);

create type lifecycle_t as enum ('active','nrnd','eol','obsolete','preview','unknown');

create table document (
  id            bigserial primary key,
  kind          text not null,                 -- datasheet | app_note | design_guide | fab_capability | standard_summary | errata | other
  title         text,
  publisher     text,
  url           text,
  local_path    text,
  sha256        char(64) unique,
  licence       text,                          -- recorded from source terms
  licence_note  text,
  fetched_at    timestamptz not null default now()
);

create table part (
  id              bigserial primary key,
  mpn             text not null,
  manufacturer_id bigint not null references manufacturer(id),
  category        text,                        -- taxonomy path 'ic/power/dc-dc/buck'
  description     text,
  package_id      bigint references package(id),
  lifecycle       lifecycle_t not null default 'unknown',
  lifecycle_src   jsonb,                       -- [{source, status, fetched_at, url}]
  rohs            boolean,
  aec_q           text,
  datasheet_id    bigint references document(id),
  parametrics     jsonb not null default '{}',
  updated_at      timestamptz not null default now(),
  unique (mpn, manufacturer_id)
);
create index on part using gin (mpn gin_trgm_ops);

create table pin (
  part_id     bigint not null references part(id) on delete cascade,
  number      text not null,
  name        text not null,
  type        text not null,                   -- power_in | power_out | input | output | bidir | passive | nc | open_collector | ...
  domain      text,                            -- supply domain / bank
  notes       text,
  provenance  jsonb not null,                  -- {doc_id, page, model, prompt_version, extracted_at}
  primary key (part_id, number)
);

create table footprint (
  id           bigserial primary key,
  part_id      bigint not null references part(id) on delete cascade,
  library_ref  text not null,                  -- 'Package_DFN_QFN:QFN-48-1EP_7x7mm_P0.5mm_EP5.6x5.6mm'
  pad_count    int not null,
  checksum     char(64),
  verified_by  text not null default 'none',   -- none | ipc_generator | vendor | human
  verified_at  timestamptz
);

create table symbol (
  id           bigserial primary key,
  part_id      bigint not null references part(id) on delete cascade,
  library_ref  text not null,
  checksum     char(64),
  verified_by  text not null default 'none',
  verified_at  timestamptz
);

create table alternate (
  part_id      bigint not null references part(id) on delete cascade,
  alt_part_id  bigint not null references part(id) on delete cascade,
  kind         text not null,                  -- drop_in | pin_compatible | functional
  notes        text,
  provenance   jsonb not null,
  primary key (part_id, alt_part_id)
);

create table offer (
  id           bigserial primary key,
  part_id      bigint not null references part(id) on delete cascade,
  distributor  text not null,
  sku          text,
  region       text not null,                  -- americas | emea | apac
  stock        int,
  moq          int,
  price_breaks jsonb,                          -- [{qty, unit_price, currency}]
  lead_days    int,
  fetched_at   timestamptz not null default now()
);
create index on offer (part_id, fetched_at desc);

create table assembler_part (
  part_id      bigint not null references part(id) on delete cascade,
  assembler    text not null,                  -- jlcpcb | pcbway | ...
  assembler_pn text not null,
  class        text,                           -- basic | extended | preferred
  stock        int,
  fetched_at   timestamptz not null default now(),
  primary key (part_id, assembler)
);

-- Document sections and chunks (store 2.2)
create table section (
  id          bigserial primary key,
  doc_id      bigint not null references document(id) on delete cascade,
  path        text not null,                   -- '8.2.3 Layout Guidelines'
  page_start  int,
  page_end    int,
  text        text not null,
  tables      jsonb not null default '[]',
  figures     text[] not null default '{}'
);
create index on section using gin (to_tsvector('english', text));

create table chunk (
  id          bigserial primary key,
  section_id  bigint not null references section(id) on delete cascade,
  text        text not null,
  token_count int,
  embedding   vector(1024)                     -- dimension set per chosen embedding model
);
create index on chunk using hnsw (embedding vector_cosine_ops);

-- Fab capabilities (store 2.5)
create table fab_capability (
  id          bigserial primary key,
  fab         text not null,
  process     text not null,
  constraints jsonb not null,                  -- min_trace_mm, min_space_mm, min_drill_mm, annular_ring_mm, ...
  stackups    jsonb not null default '[]',
  source_url  text,
  verified_at timestamptz,
  unique (fab, process)
);
