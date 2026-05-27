%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precision QPI Comparison Pipeline
%
% Purpose
% -------
% Build one strict cell correspondence table for all methods before any OVD
% comparison is performed.
%
% Ground truth:
%   ovd_1_6000nm.mat   -> 9 sub-cells, canonical indices  1:9
%   ovd_1_15000nm.mat  -> 9 sub-cells, canonical indices 10:18
%
% Methods:
%   tomo / dhm / fpm stacks can contain a larger pool, e.g. 45 sub-cells.
%   The pipeline selects the 18 cells corresponding to the two GT OVD files.
%   tie / dpc stacks contain 18 sub-cells whose import order can be wrong,
%   so they are matched to the same canonical ground-truth order.
%
% Main precision fixes:
%   1. Canonical GT order is defined once by sorted GT centroids.
%   2. Every method is solved as an exact one-to-one assignment problem.
%   3. Matching score combines local NCC, amplitude consistency, and area.
%   4. Each accepted pair is locally fine-registered before saving.
%   5. Output stacks are strictly 401 x 401 x 18 and index-compatible.
%
% Requirements:
%   Image Processing Toolbox functions used by the original workflow:
%   imresize, imrotate, normxcorr2, regionprops, padarray, bwareafilt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% -------------------- User configuration --------------------------------
cfg = struct();

% Ground-truth OVD files. The script also checks ./stacks automatically.
cfg.gt_6k_file = 'ovd_1_6000nm.mat';
cfg.gt_15k_file = 'ovd_1_15000nm.mat';

% Output folder.
cfg.output_dir = fullfile('.', 'Analysis_Results');

% Aggregate conference data file. TOMO/DHM/FPM can be stored inside this one
% .mat file instead of as separate *_comparison_data.mat files.
data_file = 'data-conference-paper.mat';
cfg.conference_data_file = fullfile('.', data_file);

% Defaults from Master_Data_Comparison_8778.m. Values in
% data-conference-paper.mat override these when present.
cfg.default_lambda_um = 0.532000000000000;
cfg.default_dx_um = 0.65;
cfg.default_polymerRI = 1.566760000000000;
cfg.default_n_immersion = 1.518431000000000;

% Strict comparison patch size.
cfg.patch_size = [401, 401];

% Segmentation and registration parameters.
cfg.gt_threshold = 0.05;
cfg.min_cell_area = 5000;
cfg.score_size = [128, 128];       % Downsampled size for fast all-pair scoring.
cfg.assignment_margin = 35;        % Full-resolution local shift search, pixels.
cfg.keep_only_gt_mask = false;     % Preserve measured cell support; OVD uses an evaluation mask.
cfg.output_mask_dilate_px = 8;     % Used only if keep_only_gt_mask is set to true.
cfg.ovd_mask_dilate_px = 8;        % Prevent strict GT masks from clipping method cells.
cfg.background_percentile = 5;     % Robust background subtraction percentile.

% DPC raw phase-map extraction. If raw maps are available, this produces a
% fresh DPC stack; otherwise the script falls back to an existing DPC stack.
cfg.run_dpc_raw_extraction_if_available = true;
cfg.dpc_base_folder = fullfile('.', 'DPC', 'results');
cfg.dpc_fov_list = [2, 6];         % 2.mat -> 6000nm, 6.mat -> 15000nm
cfg.dpc_fov_labels = {'6000', '15000'};
cfg.dpc_interactive_alignment = true;
cfg.dpc_alignment_cache = fullfile(cfg.output_dir, 'dpc_alignment_cache.mat');

% Spatial parameters for raw DPC extraction.
cfg.Mag = 10;
cfg.camera_pixel = 6.5e-6;
cfg.lambda_um = 520e-9 * 1e6;
cfg.dn = 1.55560 - 1.51041;
cfg.ROTATE_DPC_180 = true;
cfg.FLIP_RAW_DATA = true;
cfg.ROTATE_RAW_90 = true;
cfg.initial_angle_guess = 1.5;

% Candidate files for method stacks. The conference data file follows the
% convention used in Master_Data_Comparison_8778.m:
%   loaded_data.<method>.ph_stack or loaded_data.<method>.stack
% with method names such as tomo, dhmnjust, dhmwut, fpmnjust, fpmqcilab...
% If a file contains several methods, the loader first searches recursively
% for variables/fields matching the requested method name, then falls back
% to common stack variable names.
method_specs = [
    struct('name', 'tomo', 'candidate_files', {{ ...
        fullfile(cfg.output_dir, 'tomo_comparison_data.mat'), ...
        fullfile('.', 'tomo_comparison_data.mat'), ...
        fullfile('.', 'tomo.mat'), ...
        cfg.conference_data_file}})
    struct('name', 'dhmnjust', 'candidate_files', {{ ...
        fullfile(cfg.output_dir, 'dhmnjust_comparison_data.mat'), ...
        fullfile('.', 'dhmnjust_comparison_data.mat'), ...
        fullfile('.', 'dhmnjust.mat'), ...
        cfg.conference_data_file}})
    struct('name', 'dhmwut', 'candidate_files', {{ ...
        fullfile(cfg.output_dir, 'dhmwut_comparison_data.mat'), ...
        fullfile('.', 'dhmwut_comparison_data.mat'), ...
        fullfile('.', 'dhmwut.mat'), ...
        cfg.conference_data_file}})
    struct('name', 'fhpmnjust', 'candidate_files', {{ ...
        fullfile(cfg.output_dir, 'fhpmnjust_comparison_data.mat'), ...
        fullfile('.', 'fhpmnjust_comparison_data.mat'), ...
        fullfile('.', 'fhpmnjust.mat'), ...
        cfg.conference_data_file}})
    struct('name', 'fhpmwut', 'candidate_files', {{ ...
        fullfile(cfg.output_dir, 'fhpmwut_comparison_data.mat'), ...
        fullfile('.', 'fhpmwut_comparison_data.mat'), ...
        fullfile('.', 'fhpmwut.mat'), ...
        cfg.conference_data_file}})
    struct('name', 'fpmnjust', 'candidate_files', {{ ...
        fullfile(cfg.output_dir, 'fpmnjust_comparison_data.mat'), ...
        fullfile('.', 'fpmnjust_comparison_data.mat'), ...
        fullfile('.', 'fpmnjust.mat'), ...
        cfg.conference_data_file}})
    struct('name', 'fpmqcilab', 'candidate_files', {{ ...
        fullfile(cfg.output_dir, 'fpmqcilab_comparison_data.mat'), ...
        fullfile('.', 'fpmqcilab_comparison_data.mat'), ...
        fullfile('.', 'fpmqcilab.mat'), ...
        cfg.conference_data_file}})
    struct('name', 'fpmwut', 'candidate_files', {{ ...
        fullfile(cfg.output_dir, 'fpmwut_comparison_data.mat'), ...
        fullfile('.', 'fpmwut_comparison_data.mat'), ...
        fullfile('.', 'fpmwut.mat'), ...
        cfg.conference_data_file}})
    struct('name', 'tie', 'candidate_files', {{ ...
        fullfile(cfg.output_dir, 'tie_comparison_data.mat'), ...
        fullfile('.', 'tie_comparison_data.mat'), ...
        fullfile('.', 'tie.mat'), ...
        cfg.conference_data_file}})
    struct('name', 'dpc', 'candidate_files', {{ ...
        fullfile(cfg.output_dir, 'dpc_comparison_data.mat'), ...
        fullfile('.', 'dpc_comparison_data.mat'), ...
        fullfile('.', 'dpc.mat'), ...
        cfg.conference_data_file}})
];

%% -------------------- Resolve inputs ------------------------------------
if ~exist(cfg.output_dir, 'dir')
    mkdir(cfg.output_dir);
end

cfg.gt_6k_file = resolveGtPath(cfg.gt_6k_file);
cfg.gt_15k_file = resolveGtPath(cfg.gt_15k_file);

fprintf('================ Phase 1: Build canonical GT order ================\n');
gt = buildCanonicalGroundTruth(cfg);
gt.physical_params = loadGtPhysicalParams(cfg);
save(fullfile(cfg.output_dir, 'canonical_gt_401x401x18.mat'), 'gt');
fprintf('  -> Canonical GT saved: %s\n', fullfile(cfg.output_dir, 'canonical_gt_401x401x18.mat'));

%% -------------------- Optional DPC extraction ----------------------------
if cfg.run_dpc_raw_extraction_if_available && dpcRawFilesAvailable(cfg)
    fprintf('\n================ Phase 2: Extract DPC from raw phase maps ===========\n');
    dpc = struct();
    dpc.ph_stack = extractDpcStackFromRaw(cfg, gt);
    dpc.raw_note = 'Generated by precision_qpi_comparison_pipeline.m';
    save(fullfile(cfg.output_dir, 'dpc_comparison_data.mat'), 'dpc');
    fprintf('  -> Fresh DPC stack saved: %s\n', fullfile(cfg.output_dir, 'dpc_comparison_data.mat'));
else
    fprintf('\n================ Phase 2: DPC raw extraction skipped ================\n');
    fprintf('  -> Raw DPC files were not found or extraction is disabled.\n');
end

%% -------------------- Method correspondence and fine registration --------
fprintf('\n================ Phase 3: Match methods to canonical GT =============\n');
aligned = struct();
aligned.gt = gt;
method_report = struct([]);
shared_source_map = loadConferenceTruthSourceMap(cfg, gt);
shared_source_map_source = '';
if ~isempty(shared_source_map)
    shared_source_map_source = 'model_synt';
    fprintf('  -> Shared 45-to-18 source map initialized from model_synt.\n');
end

for m = 1:numel(method_specs)
    method_name = method_specs(m).name;
    [method_file, found] = firstExistingFile(method_specs(m).candidate_files);
    if ~found
        fprintf('  -> %-10s skipped: no input file found.\n', upper(method_name));
        continue;
    end

    fprintf('  -> %-10s loading: %s\n', upper(method_name), method_file);
    try
        [raw_stack, method_params] = loadMethodStack(method_file, method_name);
    catch ME
        fprintf('     skipped: %s\n', ME.message);
        continue;
    end
    raw_stack = standardizeStackSize(raw_stack, cfg.patch_size);
    method_params = fillMissingMethodParams(method_params, gt.physical_params);
    raw_stack = rescaleStackToGtPixel(raw_stack, method_params, gt.physical_params);

    if usesSharedSourceMap(method_name) && ~isempty(shared_source_map) && max(shared_source_map) <= size(raw_stack, 3)
        [ordered_stack, report] = alignStackWithSourceMap(raw_stack, gt, cfg, method_name, shared_source_map);
        report.shared_map_source = shared_source_map_source;
    else
        [ordered_stack, report] = matchAndRegisterStack(raw_stack, gt, cfg, method_name);
        if isTomoAnchorMethod(method_name)
            shared_source_map = report.selected_source_frames;
            shared_source_map_source = method_name;
            fprintf('     -> Shared source map updated from TOMO matching.\n');
        end
    end
    aligned.(method_name).ph_stack = ordered_stack;
    aligned.(method_name).source_file = method_file;
    aligned.(method_name).assignment = report.assignment;
    aligned.(method_name).score_matrix = report.score_matrix;
    aligned.(method_name).pair_scores = report.pair_scores;
    aligned.(method_name).local_shifts = report.local_shifts;
    aligned.(method_name).method_params = method_params;
    aligned.(method_name).selected_source_frames = report.selected_source_frames;

    out_struct = aligned.(method_name); %#ok<NASGU>
    out_path = fullfile(cfg.output_dir, sprintf('%s_aligned_to_gt_401x401x18.mat', method_name));
    save(out_path, 'out_struct');

    method_report = [method_report; report]; %#ok<AGROW>
    fprintf('     saved: %s | median pair score %.4f\n', out_path, median(report.pair_scores));
end

ovd_report = computeOvdReport(aligned, method_specs, cfg);
save(fullfile(cfg.output_dir, 'all_methods_aligned_to_gt_401x401x18.mat'), ...
    'aligned', 'method_report', 'ovd_report', '-v7.3');
save(fullfile(cfg.output_dir, 'ovd_comparison_report.mat'), 'ovd_report');
writeOvdCsv(ovd_report, fullfile(cfg.output_dir, 'ovd_comparison_report.csv'));
makeOvdPlots(ovd_report, cfg);

%% -------------------- Visual verification --------------------------------
fprintf('\n================ Phase 4: Verification gallery =====================\n');
makeVerificationGallery(aligned, method_specs, cfg);
fprintf('\nSUCCESS: all available method stacks are aligned to canonical GT order.\n');
fprintf('Canonical order: 1:9 = 6000nm, 10:18 = 15000nm, each 401x401.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resolved = resolveGtPath(file_name)
    if isfile(file_name)
        resolved = file_name;
        return;
    end

    stack_path = fullfile('stacks', file_name);
    if isfile(stack_path)
        resolved = stack_path;
        return;
    end

    error('Ground-truth file not found: %s', file_name);
end

function [file_name, found] = firstExistingFile(candidates)
    file_name = '';
    found = false;
    for i = 1:numel(candidates)
        if isfile(candidates{i})
            file_name = candidates{i};
            found = true;
            return;
        end
    end
end

function tf = dpcRawFilesAvailable(cfg)
    tf = true;
    for i = 1:numel(cfg.dpc_fov_list)
        mat_path = fullfile(cfg.dpc_base_folder, sprintf('%d.mat', cfg.dpc_fov_list(i)));
        tf = tf && isfile(mat_path);
    end
end

function gt = buildCanonicalGroundTruth(cfg)
    h6 = loadHeightMap(cfg.gt_6k_file);
    h15 = loadHeightMap(cfg.gt_15k_file);

    [stack6, meta6] = extractGtCells(h6, '6000nm', cfg);
    [stack15, meta15] = extractGtCells(h15, '15000nm', cfg);

    gt = struct();
    gt.stack = cat(3, stack6, stack15);
    gt.mask_stack = gt.stack > cfg.gt_threshold;
    gt.label = [repmat({'6000nm'}, 1, 9), repmat({'15000nm'}, 1, 9)];
    gt.source_file = {cfg.gt_6k_file, cfg.gt_15k_file};
    gt.patch_size = cfg.patch_size;
    gt.meta = [meta6, meta15];

    for k = 1:18
        gt.index_table(k).canonical_index = k; %#ok<AGROW>
        gt.index_table(k).ovd = gt.label{k}; %#ok<AGROW>
        gt.index_table(k).source_cell_index = gt.meta(k).source_cell_index; %#ok<AGROW>
        gt.index_table(k).centroid = gt.meta(k).centroid; %#ok<AGROW>
    end
end

function params = loadGtPhysicalParams(cfg)
    params = struct();
    params.lambda = cfg.default_lambda_um;
    params.dx = cfg.default_dx_um;
    params.polymerRI = cfg.default_polymerRI;
    params.n_immersion = cfg.default_n_immersion;

    if isfile(cfg.conference_data_file)
        s = load(cfg.conference_data_file);
        if isfield(s, 'model_synt')
            params = mergeMethodParams(params, extractMethodParams(s.model_synt));
        end
    end
    params.delta_n = params.polymerRI - params.n_immersion;
end

function source_map = loadConferenceTruthSourceMap(cfg, gt)
    source_map = [];
    if ~isfile(cfg.conference_data_file)
        return;
    end

    try
        s = load(cfg.conference_data_file);
        if ~isfield(s, 'model_synt')
            return;
        end
        model_stack = unwrapStackVariable(s.model_synt);
        if isempty(model_stack) || size(model_stack, 3) < 18
            return;
        end
        model_stack = standardizeStackSize(model_stack, cfg.patch_size);
        source_patches = buildCandidatePatchStack(model_stack, cfg);
        score_matrix = computeScoreMatrix(source_patches, gt, cfg);
        source_map = exactMaxAssignment(score_matrix);
    catch ME
        fprintf('  -> model_synt source map skipped: %s\n', ME.message);
        source_map = [];
    end
end

function tf = isTomoAnchorMethod(method_name)
    tf = strcmpi(method_name, 'tomo');
end

function tf = usesSharedSourceMap(method_name)
    lower_name = lower(method_name);
    tf = strcmp(lower_name, 'tomo') || ~isempty(strfind(lower_name, 'dhm')) || ...
        ~isempty(strfind(lower_name, 'fpm'));
end

function height_map = loadHeightMap(mat_path)
    s = load(mat_path);
    if isfield(s, 'height_map')
        height_map = double(s.height_map);
        return;
    end

    names = fieldnames(s);
    for i = 1:numel(names)
        v = s.(names{i});
        if isnumeric(v) && ndims(v) == 2
            height_map = double(v);
            return;
        end
    end
    error('No 2-D height_map-like variable found in %s', mat_path);
end

function [stack, meta] = extractGtCells(height_map, label, cfg)
    bw = height_map > cfg.gt_threshold;
    bw = bwareafilt(bw, 9);
    props = regionprops(bw, height_map, 'Centroid', 'Area', 'BoundingBox', 'MaxIntensity');
    props = props([props.Area] >= cfg.min_cell_area);
    if numel(props) < 9
        error('%s has only %d detected cells; expected 9.', label, numel(props));
    end

    [~, area_order] = sort([props.Area], 'descend');
    props = props(area_order(1:9));
    order = canonicalGridOrder(cat(1, props.Centroid));
    props = props(order);

    H = cfg.patch_size(1);
    W = cfg.patch_size(2);
    stack = zeros(H, W, 9);
    meta = repmat(struct('label', '', 'source_cell_index', [], 'centroid', [], ...
        'area', [], 'bounding_box', [], 'max_intensity', []), 1, 9);

    for i = 1:9
        center_xy = props(i).Centroid;
        patch = cropCenteredWithPadding(height_map, center_xy(2), center_xy(1), H, W, 0);
        mask = patch > cfg.gt_threshold;
        if any(mask(:))
            mask = bwareafilt(mask, 1);
            patch(~mask) = 0;
            patch = robustZeroPatch(patch, mask, cfg.background_percentile);
        end

        stack(:, :, i) = patch;
        meta(i).label = label;
        meta(i).source_cell_index = i;
        meta(i).centroid = center_xy;
        meta(i).area = props(i).Area;
        meta(i).bounding_box = props(i).BoundingBox;
        meta(i).max_intensity = props(i).MaxIntensity;
    end
end

function order = canonicalGridOrder(centroids)
    % Sort a 3x3 cell array top-to-bottom, then left-to-right in each row.
    if size(centroids, 1) ~= 9
        [~, order] = sortrows(centroids, [2, 1]);
        return;
    end

    [~, y_order] = sort(centroids(:, 2), 'ascend');
    order = zeros(9, 1);
    out = 1;
    for row = 1:3
        idx = y_order((row - 1) * 3 + (1:3));
        [~, x_local] = sort(centroids(idx, 1), 'ascend');
        order(out:out + 2) = idx(x_local);
        out = out + 3;
    end
end

function stack = extractDpcStackFromRaw(cfg, gt)
    dpc_stack = zeros(cfg.patch_size(1), cfg.patch_size(2), 18);
    cache = loadDpcCache(cfg);

    for idx = 1:numel(cfg.dpc_fov_list)
        fov = cfg.dpc_fov_list(idx);
        mat_path = fullfile(cfg.dpc_base_folder, sprintf('%d.mat', fov));
        fprintf('  -> Processing DPC raw phase map: %s\n', mat_path);

        dpc_raw = loadFirstNumericMatrix(mat_path);
        if size(dpc_raw, 3) == 3
            dpc_raw = mean(dpc_raw, 3);
        end
        dpc_raw = double(dpc_raw);
        if cfg.ROTATE_DPC_180
            dpc_raw = rot90(dpc_raw, 2);
        end

        if idx == 1
            gt_file = cfg.gt_6k_file;
            target_indices = 1:9;
        else
            gt_file = cfg.gt_15k_file;
            target_indices = 10:18;
        end

        gt_data = load(gt_file);
        design_p = cfg.dn * double(gt_data.height_map) .* (2 * pi / cfg.lambda_um);
        if cfg.FLIP_RAW_DATA
            design_p = fliplr(design_p);
        end
        if cfg.ROTATE_RAW_90
            design_p = rot90(design_p, -1);
        end

        meas_sampling_um = (cfg.camera_pixel / cfg.Mag) * 1e6;
        if isfield(gt_data, 'sampling')
            raw_scale = gt_data.sampling(1) / meas_sampling_um;
            if raw_scale < 0.5 || raw_scale > 2.0
                raw_scale = 1.0;
            end
        else
            raw_scale = 1.0;
        end

        cache_key = sprintf('fov_%d', fov);
        if isfield(cache, cache_key)
            params = cache.(cache_key);
        else
            params = estimateDpcAlignment(dpc_raw, design_p, raw_scale, fov, cfg);
            cache.(cache_key) = params; %#ok<AGROW>
            saveDpcCache(cfg, cache);
        end

        ph_corrected = removeQuadraticBackground(dpc_raw, params, design_p, cfg);
        ph_native = inverseMapDpcToNative(ph_corrected, params, design_p, cfg);

        group_stack = extractNativeCellsWithGtOrder(ph_native, gt, target_indices, cfg);
        dpc_stack(:, :, target_indices) = group_stack;
    end

    [dpc_stack, ~] = matchAndRegisterStack(dpc_stack, gt, cfg, 'dpc_raw_verified');
end

function cache = loadDpcCache(cfg)
    cache = struct();
    if isfile(cfg.dpc_alignment_cache)
        tmp = load(cfg.dpc_alignment_cache);
        if isfield(tmp, 'cache')
            cache = tmp.cache;
        end
    end
end

function saveDpcCache(cfg, cache)
    if ~exist(cfg.output_dir, 'dir')
        mkdir(cfg.output_dir);
    end
    save(cfg.dpc_alignment_cache, 'cache');
end

function params = estimateDpcAlignment(dpc_full_phase, design_p, raw_scale, fov, cfg)
    design_unrotated = imresize(design_p, raw_scale, 'bilinear');

    if cfg.dpc_interactive_alignment
        f_crop = figure('Name', sprintf('FOV %d - draw rough ROI', fov), ...
            'Units', 'normalized', 'OuterPosition', [0.1 0.1 0.8 0.8]);
        imagesc(dpc_full_phase); axis image; colormap gray;
        title(sprintf('FOV %d: draw ROI around the 9-cell array, then double-click inside it.', fov), ...
            'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');
        roi_obj = drawrectangle(gca, 'Color', 'r', 'Label', 'rough ROI');
        wait(roi_obj);
        roi = round(roi_obj.Position);
        close(f_crop);
    else
        roi = [1, 1, size(dpc_full_phase, 2), size(dpc_full_phase, 1)];
    end

    ph_crop = imcrop(dpc_full_phase, roi);
    mask_gt = design_unrotated > cfg.gt_threshold;
    props = regionprops(mask_gt, 'BoundingBox');
    if ~isempty(props)
        bbox_all = cat(1, props.BoundingBox);
        gt_w = max(bbox_all(:, 1) + bbox_all(:, 3)) - min(bbox_all(:, 1));
        gt_h = max(bbox_all(:, 2) + bbox_all(:, 4)) - min(bbox_all(:, 2));
        rough_scale = mean([roi(3) / gt_w, roi(4) / gt_h]);
    else
        rough_scale = 1.0;
    end

    angle = cfg.initial_angle_guess;
    scales_to_test = rough_scale * linspace(0.75, 1.25, 31);
    best_ncc = -inf;
    best_scale = rough_scale;
    best_shift = [0, 0];

    fprintf('     coarse scale search around %.3f ...\n', rough_scale);
    target_small = normalizeForNcc(imresize(ph_crop, 0.5, 'bilinear'));
    for s_idx = 1:numel(scales_to_test)
        design_test = imrotate(imresize(design_unrotated, scales_to_test(s_idx), 'bilinear'), ...
            angle, 'bilinear', 'crop');
        design_small = normalizeForNcc(imresize(design_test, 0.5, 'bilinear'));
        target_padded = padTargetIfNeeded(target_small, size(design_small));
        c = normxcorr2(design_small, target_padded);
        [max_c, imax] = max(c(:));
        if max_c > best_ncc
            best_ncc = max_c;
            best_scale = scales_to_test(s_idx);
            [ypeak, xpeak] = ind2sub(size(c), imax);
            best_shift = 2 * [ypeak - size(design_small, 1), xpeak - size(design_small, 2)];
        end
    end

    curr_angle = angle;
    curr_scale = best_scale;
    total_shift = best_shift;
    curr_flip = false;

    if cfg.dpc_interactive_alignment
        fig_reg = figure('Name', sprintf('FOV %d alignment control', fov), ...
            'Units', 'normalized', 'OuterPosition', [0.1 0.1 0.8 0.8]);
        done = false;
        while ~done && ishandle(fig_reg)
            if curr_flip
                design_t = fliplr(design_unrotated);
            else
                design_t = design_unrotated;
            end
            design_r = imrotate(imresize(design_t, curr_scale, 'bilinear'), ...
                curr_angle, 'bilinear', 'crop');

            RGB = makeOverlay(ph_crop, design_r, total_shift);
            imshow(RGB);
            title(sprintf(['FOV %d: ENTER confirm | angle %.2f | scale %.4f | shift [%d,%d]\n' ...
                'W/A/S/D move | Q/E rotate | Z/X scale | F flip'], ...
                fov, curr_angle, curr_scale, round(total_shift(1)), round(total_shift(2))), ...
                'FontSize', 13, 'FontWeight', 'bold');

            if waitforbuttonpress == 1
                kp = double(get(fig_reg, 'CurrentCharacter'));
                if kp == 13
                    done = true;
                elseif kp == 97 || kp == 28
                    total_shift(2) = total_shift(2) - 1;
                elseif kp == 100 || kp == 29
                    total_shift(2) = total_shift(2) + 1;
                elseif kp == 119 || kp == 30
                    total_shift(1) = total_shift(1) - 1;
                elseif kp == 115 || kp == 31
                    total_shift(1) = total_shift(1) + 1;
                elseif kp == 113
                    curr_angle = curr_angle + 0.1;
                elseif kp == 101
                    curr_angle = curr_angle - 0.1;
                elseif kp == 102
                    curr_flip = ~curr_flip;
                elseif kp == 122 || kp == 45
                    curr_scale = curr_scale - 0.005;
                elseif kp == 120 || kp == 61 || kp == 43
                    curr_scale = curr_scale + 0.005;
                end
            end
        end
        if ishandle(fig_reg)
            close(fig_reg);
        end
    end

    params = struct();
    params.roi = roi;
    params.global_shift = [roi(2) - 1 + total_shift(1), roi(1) - 1 + total_shift(2)];
    params.angle = curr_angle;
    params.scale = curr_scale;
    params.flip = curr_flip;
    params.raw_scale = raw_scale;
end

function img = padTargetIfNeeded(img, template_size)
    pad_y = max(0, template_size(1) - size(img, 1));
    pad_x = max(0, template_size(2) - size(img, 2));
    if pad_y > 0 || pad_x > 0
        img = padarray(img, [pad_y, pad_x], mean(img(:)), 'post');
    end
end

function RGB = makeOverlay(measured, reference, shift_yx)
    [h_m, w_m] = size(measured);
    [h_r, w_r] = size(reference);
    h_c = max(h_m, h_r) + 200;
    w_c = max(w_m, w_r) + 200;

    m_p = padarray(measured, [h_c - h_m, w_c - w_m], 'replicate', 'post');
    r_p = padarray(reference, [h_c - h_r, w_c - w_r], 0, 'post');

    nm = normalize01(m_p);
    nr = normalize01(r_p);
    rs = imtranslate(nr, [round(shift_yx(2)), round(shift_yx(1))], 'OutputView', 'same');

    RGB = zeros(size(nm, 1), size(nm, 2), 3);
    RGB(:, :, 2) = nm;
    RGB(:, :, 1) = rs;
    RGB(:, :, 3) = rs;
end

function ph_corrected = removeQuadraticBackground(dpc_full_phase, params, design_p, cfg)
    [H_f, W_f] = size(dpc_full_phase);
    design_unrotated = imresize(design_p, params.raw_scale, 'bilinear');
    if params.flip
        design_unrotated = fliplr(design_unrotated);
    end
    final_design = imrotate(imresize(design_unrotated, params.scale, 'bilinear'), ...
        params.angle, 'bilinear', 'crop');

    mask_bgr = true(H_f, W_f);
    y0 = round(params.global_shift(1));
    x0 = round(params.global_shift(2));
    r_start_y = max(1, y0 + 1);
    r_end_y = min(H_f, y0 + size(final_design, 1));
    r_start_x = max(1, x0 + 1);
    r_end_x = min(W_f, x0 + size(final_design, 2));
    mask_bgr(r_start_y:r_end_y, r_start_x:r_end_x) = false;

    [Y_grid, X_grid] = ndgrid(1:H_f, 1:W_f);
    A_fit = [X_grid(mask_bgr).^2, Y_grid(mask_bgr).^2, X_grid(mask_bgr).*Y_grid(mask_bgr), ...
        X_grid(mask_bgr), Y_grid(mask_bgr), ones(sum(mask_bgr(:)), 1)];
    coeffs = A_fit \ dpc_full_phase(mask_bgr);
    bg = coeffs(1) * X_grid.^2 + coeffs(2) * Y_grid.^2 + coeffs(3) * X_grid .* Y_grid + ...
        coeffs(4) * X_grid + coeffs(5) * Y_grid + coeffs(6);
    ph_corrected = dpc_full_phase - bg;
end

function ph_native = inverseMapDpcToNative(ph_corrected, params, design_p, cfg)
    [H_f, W_f] = size(ph_corrected);
    design_unrotated = imresize(design_p, params.raw_scale, 'bilinear');
    if params.flip
        design_t = fliplr(design_unrotated);
    else
        design_t = design_unrotated;
    end
    design_r = imrotate(imresize(design_t, params.scale, 'bilinear'), params.angle, 'bilinear', 'crop');
    H_u = size(design_r, 1);
    W_u = size(design_r, 2);

    ph_box = zeros(H_u, W_u);
    y0 = round(params.global_shift(1));
    x0 = round(params.global_shift(2));
    y_d_s = max(1, 1 + y0);
    y_d_e = min(H_f, H_u + y0);
    x_d_s = max(1, 1 + x0);
    x_d_e = min(W_f, W_u + x0);
    y_b_s = max(1, 1 - y0);
    y_b_e = min(H_u, H_f - y0);
    x_b_s = max(1, 1 - x0);
    x_b_e = min(W_u, W_f - x0);
    ph_box(y_b_s:y_b_e, x_b_s:x_b_e) = ph_corrected(y_d_s:y_d_e, x_d_s:x_d_e);

    ph_unrotated = imrotate(ph_box, -params.angle, 'bilinear', 'crop');
    if params.flip
        ph_unrotated = fliplr(ph_unrotated);
    end

    ph_design_p = imresize(ph_unrotated, size(design_p), 'bilinear');
    ph_native = ph_design_p;
    if cfg.ROTATE_RAW_90
        ph_native = rot90(ph_native, 1);
    end
    if cfg.FLIP_RAW_DATA
        ph_native = fliplr(ph_native);
    end
end

function stack = extractNativeCellsWithGtOrder(ph_native, gt, target_indices, cfg)
    H = cfg.patch_size(1);
    W = cfg.patch_size(2);
    stack = zeros(H, W, numel(target_indices));
    pad_amt = 1500;
    ph_padded = padarray(ph_native, [pad_amt, pad_amt], 0, 'both');

    for j = 1:numel(target_indices)
        k = target_indices(j);
        cx = round(gt.meta(k).centroid(1)) + pad_amt;
        cy = round(gt.meta(k).centroid(2)) + pad_amt;
        raw_patch = cropCenteredWithPadding(ph_padded, cy, cx, H, W, 0);
        target = gt.stack(:, :, k);
        mask = gt.mask_stack(:, :, k);
        [aligned_patch, ~] = localRegisterPatch(raw_patch, target, mask, cfg.assignment_margin);
        aligned_patch = cleanWithMask(aligned_patch, mask, cfg);
        stack(:, :, j) = aligned_patch;
    end
end

function data = loadFirstNumericMatrix(mat_path)
    s = load(mat_path);
    names = fieldnames(s);
    for i = 1:numel(names)
        v = s.(names{i});
        if isnumeric(v) && ndims(v) >= 2
            data = double(v);
            return;
        end
    end
    error('No numeric matrix found in %s', mat_path);
end

function [stack, method_params] = loadMethodStack(mat_path, method_name)
    s = load(mat_path);
    aliases = methodAliases(method_name);
    method_params = defaultMethodParams(method_name, mat_path);

    % Aggregate files can contain several methods. Prefer a stack whose
    % variable path explicitly names the requested method.
    [stack, source_path, method_params] = findNamedStack(s, aliases, '', method_params);
    if ~isempty(stack)
        fprintf('     method stack found at variable path: %s\n', source_path);
        return;
    end

    if strcmpi(fileNameOnly(mat_path), 'data-conference-paper.mat')
        error('Method "%s" was not found as a named field/path in %s.', method_name, mat_path);
    end

    preferred = {'ph_stack', 'stack', 'data', 'phase_stack', method_name};
    for i = 1:numel(preferred)
        if isfield(s, preferred{i})
            stack = unwrapStackVariable(s.(preferred{i}));
            if ~isempty(stack)
                method_params = mergeMethodParams(method_params, extractMethodParams(s.(preferred{i})));
                method_params.stack_path = preferred{i};
                return;
            end
        end
    end

    names = fieldnames(s);
    for i = 1:numel(names)
        stack = unwrapStackVariable(s.(names{i}));
        if ~isempty(stack)
            method_params = mergeMethodParams(method_params, extractMethodParams(s.(names{i})));
            method_params.stack_path = names{i};
            return;
        end
    end

    error('No usable 18-layer numeric stack found in %s', mat_path);
end

function params = defaultMethodParams(method_name, mat_path)
    params = struct();
    params.method = method_name;
    params.source_file = mat_path;
    params.stack_path = '';
    params.lambda = NaN;
    params.dx = NaN;
    params.polymerRI = NaN;
    params.n_immersion = NaN;
    params.delta_n = NaN;
end

function params = fillMissingMethodParams(params, gt_params)
    names = {'lambda', 'dx', 'polymerRI', 'n_immersion'};
    for i = 1:numel(names)
        if ~isfield(params, names{i}) || ~isfinite(params.(names{i}))
            params.(names{i}) = gt_params.(names{i});
        end
    end
    params.delta_n = params.polymerRI - params.n_immersion;
end

function name = fileNameOnly(path_name)
    [~, stem, ext] = fileparts(path_name);
    name = [stem, ext];
end

function aliases = methodAliases(method_name)
    switch lower(method_name)
        case 'tomo'
            aliases = {'tomo', 'tomography'};
        case 'dhm'
            aliases = {'dhm', 'digitalholography', 'digital_holography'};
        case 'fpm'
            aliases = {'fpm', 'fourierptychography', 'fourier_ptychography'};
        case 'tie'
            aliases = {'tie'};
        case 'dpc'
            aliases = {'dpc'};
        otherwise
            aliases = {lower(method_name)};
    end
end

function [stack, source_path, method_params] = findNamedStack(v, aliases, current_path, inherited_params)
    stack = [];
    source_path = '';
    method_params = inherited_params;

    if isNamedStackCandidate(v, aliases, current_path)
        stack = unwrapStackVariable(v);
        if ~isempty(stack)
            source_path = current_path;
            method_params = mergeMethodParams(inherited_params, extractMethodParams(v));
            method_params.stack_path = current_path;
            return;
        end
    end

    if isstruct(v)
        names = fieldnames(v);

        % First inspect fields whose names directly match the method. This
        % prevents an aggregate file from returning another method's stack.
        for pass = 1:2
            for i = 1:numel(names)
                field_path = appendPath(current_path, names{i});
                field_matches = nameMatchesMethod(names{i}, aliases);
                if (pass == 1 && ~field_matches) || (pass == 2 && field_matches)
                    continue;
                end

                next_params = inherited_params;
                if field_matches
                    next_params = mergeMethodParams(inherited_params, extractMethodParams(v.(names{i})));
                end
                [stack, source_path, method_params] = findNamedStack(v.(names{i}), aliases, field_path, next_params);
                if ~isempty(stack)
                    return;
                end
            end
        end
    elseif iscell(v)
        for i = 1:numel(v)
            item_path = sprintf('%s{%d}', current_path, i);
            [stack, source_path, method_params] = findNamedStack(v{i}, aliases, item_path, inherited_params);
            if ~isempty(stack)
                return;
            end
        end
    end
end

function tf = isNamedStackCandidate(v, aliases, current_path)
    tf = ~isempty(current_path) && nameMatchesMethod(current_path, aliases) && ...
        ((isnumeric(v) && ndims(v) == 3 && hasEnoughStackLayers(v)) || isstruct(v) || iscell(v));
end

function tf = nameMatchesMethod(name, aliases)
    normalized_name = normalizeName(name);
    tf = false;
    for i = 1:numel(aliases)
        if ~isempty(strfind(normalized_name, normalizeName(aliases{i}))) %#ok<STREMP>
            tf = true;
            return;
        end
    end
end

function out = normalizeName(in)
    out = lower(regexprep(char(in), '[^a-zA-Z0-9]', ''));
end

function out = appendPath(base_path, field_name)
    if isempty(base_path)
        out = field_name;
    else
        out = sprintf('%s.%s', base_path, field_name);
    end
end

function params = extractMethodParams(v)
    params = struct();
    if ~isstruct(v)
        return;
    end

    params.lambda = readNumericScalarField(v, {'lambda', 'lambda_um', 'wavelength'});
    params.dx = readNumericScalarField(v, {'dx', 'pixel_size', 'pixelsize', 'sampling'});
    params.polymerRI = readNumericScalarField(v, {'polymerRI', 'polymer_ri', 'n_polymer'});
    params.n_immersion = readNumericScalarField(v, {'n_immersion', 'nImmersion', 'immersionRI', 'n_medium'});
    if isfield(params, 'polymerRI') && isfield(params, 'n_immersion') && ...
            isfinite(params.polymerRI) && isfinite(params.n_immersion)
        params.delta_n = params.polymerRI - params.n_immersion;
    end
end

function value = readNumericScalarField(s, names)
    value = NaN;
    for i = 1:numel(names)
        if isfield(s, names{i})
            candidate = s.(names{i});
            if isnumeric(candidate) && ~isempty(candidate)
                value = double(candidate(1));
                return;
            end
        end
    end
end

function out = mergeMethodParams(base_params, new_params)
    out = base_params;
    names = fieldnames(new_params);
    for i = 1:numel(names)
        value = new_params.(names{i});
        if isnumeric(value)
            if isscalar(value) && isfinite(value)
                out.(names{i}) = value;
            end
        elseif ~isempty(value)
            out.(names{i}) = value;
        end
    end
    if isfield(out, 'polymerRI') && isfield(out, 'n_immersion') && ...
            isfinite(out.polymerRI) && isfinite(out.n_immersion)
        out.delta_n = out.polymerRI - out.n_immersion;
    end
end

function stack = unwrapStackVariable(v)
    stack = [];
    if isnumeric(v) && ndims(v) == 3
        if hasEnoughStackLayers(v)
            stack = double(v);
            stack = moveLayerDimToThird(stack);
        end
        return;
    end

    if isstruct(v)
        names = fieldnames(v);
        preferred = {'ph_stack', 'stack', 'data', 'phase_stack'};
        for p = 1:numel(preferred)
            if isfield(v, preferred{p})
                stack = unwrapStackVariable(v.(preferred{p}));
                if ~isempty(stack)
                    return;
                end
            end
        end
        for i = 1:numel(names)
            stack = unwrapStackVariable(v.(names{i}));
            if ~isempty(stack)
                return;
            end
        end
    end

    if iscell(v)
        for i = 1:numel(v)
            stack = unwrapStackVariable(v{i});
            if ~isempty(stack)
                return;
            end
        end
    end
end

function stack = moveLayerDimToThird(stack)
    dims = size(stack);
    if numel(dims) ~= 3
        error('Stack must be 3-D.');
    end

    layer_dim = findLayerDim(dims);
    if layer_dim == 3
        return;
    elseif layer_dim == 1
        stack = permute(stack, [2, 3, 1]);
    elseif layer_dim == 2
        stack = permute(stack, [1, 3, 2]);
    else
        error('Cannot find a valid cell-frame dimension with at least 18 layers.');
    end
end

function tf = hasEnoughStackLayers(stack)
    dims = size(stack);
    tf = numel(dims) == 3 && findLayerDim(dims) > 0;
end

function layer_dim = findLayerDim(dims)
    candidates = find(dims >= 18 & dims <= 200);
    if isempty(candidates)
        layer_dim = 0;
        return;
    end
    [~, best_idx] = min(dims(candidates));
    layer_dim = candidates(best_idx);
end

function stack = standardizeStackSize(stack, patch_size)
    stack = moveLayerDimToThird(stack);
    if size(stack, 3) < 18
        error('Expected at least 18 sub-cell layers, got %d.', size(stack, 3));
    end
end

function stack = rescaleStackToGtPixel(stack, method_params, gt_params)
    if ~isfield(method_params, 'dx') || ~isfinite(method_params.dx) || ...
            ~isfield(gt_params, 'dx') || ~isfinite(gt_params.dx) || gt_params.dx <= 0
        return;
    end

    scale_factor = method_params.dx / gt_params.dx;
    if ~isfinite(scale_factor) || scale_factor <= 0 || abs(scale_factor - 1) < 1e-3
        return;
    end

    first = imresize(double(stack(:, :, 1)), scale_factor, 'bilinear');
    out = zeros(size(first, 1), size(first, 2), size(stack, 3));
    out(:, :, 1) = first;
    for k = 2:size(stack, 3)
        out(:, :, k) = imresize(double(stack(:, :, k)), scale_factor, 'bilinear');
    end
    stack = out;
end

function [ordered_stack, report] = matchAndRegisterStack(raw_stack, gt, cfg, method_name)
    raw_stack = standardizeStackSize(raw_stack, cfg.patch_size);
    source_patches = buildCandidatePatchStack(raw_stack, cfg);
    score_matrix = computeScoreMatrix(source_patches, gt, cfg);
    source_for_target = exactMaxAssignment(score_matrix);

    [ordered_stack, local_shifts, pair_scores] = alignSelectedSourceFrames( ...
        source_patches, source_for_target, score_matrix, gt, cfg);

    report = buildAssignmentReport(method_name, source_for_target, score_matrix, pair_scores, local_shifts);
    printAssignmentReport(method_name, report, gt);
end

function [ordered_stack, report] = alignStackWithSourceMap(raw_stack, gt, cfg, method_name, source_for_target)
    raw_stack = standardizeStackSize(raw_stack, cfg.patch_size);
    source_patches = buildCandidatePatchStack(raw_stack, cfg);
    if max(source_for_target) > size(source_patches, 3)
        error('Shared source map references frame %d, but %s has only %d frames.', ...
            max(source_for_target), method_name, size(source_patches, 3));
    end

    score_matrix = computeScoreMatrix(source_patches, gt, cfg);
    [ordered_stack, local_shifts, pair_scores] = alignSelectedSourceFrames( ...
        source_patches, source_for_target, score_matrix, gt, cfg);

    report = buildAssignmentReport(method_name, source_for_target, score_matrix, pair_scores, local_shifts);
    report.used_shared_source_map = true;
    printAssignmentReport(method_name, report, gt);
end

function [ordered_stack, local_shifts, pair_scores] = alignSelectedSourceFrames(source_patches, source_for_target, score_matrix, gt, cfg)
    n_target = size(gt.stack, 3);
    ordered_stack = zeros(size(gt.stack));
    local_shifts = zeros(n_target, 2);
    pair_scores = zeros(n_target, 1);

    for dst = 1:n_target
        src = source_for_target(dst);
        mask = expandedMask(gt.mask_stack(:, :, dst), cfg.output_mask_dilate_px);
        [aligned_patch, shift_yx] = localRegisterPatch(source_patches(:, :, src), gt.stack(:, :, dst), ...
            mask, cfg.assignment_margin);
        ordered_stack(:, :, dst) = cleanWithMask(aligned_patch, mask, cfg);
        local_shifts(dst, :) = shift_yx;
        pair_scores(dst) = score_matrix(src, dst);
    end
end

function report = buildAssignmentReport(method_name, source_for_target, score_matrix, pair_scores, local_shifts)
    report = struct();
    report.method = method_name;
    report.assignment = [source_for_target(:), (1:numel(source_for_target))'];
    report.score_matrix = score_matrix;
    report.pair_scores = pair_scores;
    report.local_shifts = local_shifts;
    report.selected_source_frames = source_for_target(:)';
    report.used_shared_source_map = false;
    report.shared_map_source = '';
end

function source_patches = buildCandidatePatchStack(raw_stack, cfg)
    H = cfg.patch_size(1);
    W = cfg.patch_size(2);
    n_source = size(raw_stack, 3);
    source_patches = zeros(H, W, n_source);
    for k = 1:n_source
        source_patches(:, :, k) = extractCandidatePatch(raw_stack(:, :, k), cfg.patch_size);
    end
end

function patch = extractCandidatePatch(img, patch_size)
    H = patch_size(1);
    W = patch_size(2);
    img = double(img);
    img(~isfinite(img)) = 0;

    if size(img, 1) == H && size(img, 2) == W
        patch = img;
        return;
    end

    [cy, cx] = estimateCellCenter(img);
    patch = cropCenteredWithPadding(img, cy, cx, H, W, median(img(:)));
end

function [cy, cx] = estimateCellCenter(img)
    bg = median(img(:));
    signal = abs(img - bg);
    finite_signal = signal(isfinite(signal));
    if isempty(finite_signal) || max(finite_signal) <= 0
        cy = (size(img, 1) + 1) / 2;
        cx = (size(img, 2) + 1) / 2;
        return;
    end

    thresh = max(prctile(finite_signal, 92) * 0.35, max(finite_signal) * 0.08);
    bw = signal > thresh;
    if any(bw(:))
        bw = bwareafilt(bw, 1);
        props = regionprops(bw, 'Centroid', 'Area');
        if ~isempty(props)
            [~, idx] = max([props.Area]);
            cx = props(idx).Centroid(1);
            cy = props(idx).Centroid(2);
            return;
        end
    end

    [~, imax] = max(signal(:));
    [cy, cx] = ind2sub(size(signal), imax);
end

function score_matrix = computeScoreMatrix(raw_stack, gt, cfg)
    n_source = size(raw_stack, 3);
    n_target = size(gt.stack, 3);
    score_matrix = zeros(n_source, n_target);

    source_desc = patchDescriptors(raw_stack);
    target_desc = patchDescriptors(gt.stack);
    amp_scale = safeMedianPositive(source_desc.amp) / max(safeMedianPositive(target_desc.amp), eps);
    area_scale = safeMedianPositive(source_desc.area) / max(safeMedianPositive(target_desc.area), eps);

    source_feat = cell(n_source, 1);
    source_shape = cell(n_source, 1);
    target_feat = cell(n_target, 1);
    target_shape = cell(n_target, 1);
    for i = 1:n_source
        source_feat{i} = makeScoreFeature(raw_stack(:, :, i), cfg.score_size);
        source_shape{i} = makeShapeFeature(raw_stack(:, :, i), cfg.score_size);
    end
    for i = 1:n_target
        target_feat{i} = makeScoreFeature(gt.stack(:, :, i), cfg.score_size);
        target_shape{i} = makeShapeFeature(gt.stack(:, :, i), cfg.score_size);
    end

    for src = 1:n_source
        for dst = 1:n_target
            ncc_score = localNccScore(source_feat{src}, target_feat{dst}, 8);
            shape_score = localNccScore(source_shape{src}, target_shape{dst}, 8);

            predicted_amp = target_desc.amp(dst) * amp_scale;
            amp_score = exp(-abs(log((source_desc.amp(src) + eps) / (predicted_amp + eps))) / 0.55);

            predicted_area = target_desc.area(dst) * area_scale;
            area_score = exp(-abs(log((source_desc.area(src) + eps) / (predicted_area + eps))) / 0.80);

            score_matrix(src, dst) = 0.42 * ncc_score + 0.40 * shape_score + ...
                0.12 * amp_score + 0.06 * area_score;
        end
    end
end

function desc = patchDescriptors(stack)
    n = size(stack, 3);
    desc.amp = zeros(n, 1);
    desc.area = zeros(n, 1);
    for k = 1:n
        p = double(stack(:, :, k));
        p(~isfinite(p)) = 0;
        nz = abs(p) > max(1e-9, 0.03 * max(abs(p(:))));
        desc.area(k) = sum(nz(:));
        vals = p(nz);
        if isempty(vals)
            vals = p(:);
        end
        vals = vals(isfinite(vals));
        if isempty(vals)
            desc.amp(k) = 0;
        else
            desc.amp(k) = prctile(abs(vals), 95);
        end
    end
end

function feat = makeScoreFeature(patch, score_size)
    patch = double(patch);
    patch(~isfinite(patch)) = 0;
    patch = robustZeroPatch(patch, abs(patch) > max(1e-9, 0.03 * max(abs(patch(:)))), 5);
    if max(abs(patch(:))) > 0
        patch = patch / max(abs(patch(:)));
    end
    feat = imresize(patch, score_size, 'bilinear');
    feat = normalizeForNcc(feat);
end

function feat = makeShapeFeature(patch, score_size)
    patch = double(patch);
    mask = adaptiveCellMask(patch);
    feat = imresize(double(mask), score_size, 'bilinear');
    feat = normalizeForNcc(feat);
end

function mask = adaptiveCellMask(patch)
    patch = double(patch);
    patch(~isfinite(patch)) = 0;
    bg = median(patch(:));
    signal = abs(patch - bg);
    max_signal = max(signal(:));
    if max_signal <= 0
        mask = false(size(patch));
        return;
    end
    thresh = max(prctile(signal(:), 88) * 0.40, max_signal * 0.06);
    mask = signal > thresh;
    if any(mask(:))
        mask = bwareafilt(mask, 1);
    end
end

function score = localNccScore(source_feat, target_feat, margin)
    search = padarray(source_feat, [margin, margin], 0, 'both');
    c = normxcorr2(target_feat, search);
    score = max(c(:));
    if ~isfinite(score)
        score = -1;
    end
end

function assignment = exactMaxAssignment(score_matrix)
    % Rectangular exact assignment. Rows are source frames; columns are the
    % 18 canonical GT cells. Source frames may be skipped, so a 45-frame
    % TOMO/DHM/FPM stack is reduced to the best 18 one-to-one matches.
    [n_source, n_target] = size(score_matrix);
    if n_source < n_target
        error('Need at least %d source frames, got %d.', n_target, n_source);
    end

    if exist('matchpairs', 'file') == 2
        try
            pairs = matchpairs(-score_matrix, 1e6);
            if size(pairs, 1) >= n_target && numel(unique(pairs(:, 2))) == n_target
                assignment = zeros(n_target, 1);
                for i = 1:size(pairs, 1)
                    assignment(pairs(i, 2)) = pairs(i, 1);
                end
                if all(assignment > 0)
                    return;
                end
            end
        catch
            % Fall back to the DP implementation below.
        end
    end

    if n_target > 22
        error('Exact DP assignment is intended for <=22 target cells; got %d.', n_target);
    end

    num_states = 2^n_target;
    dp = -inf(num_states, 1);
    parent_col = zeros(num_states, n_source, 'uint16');
    parent_state = zeros(num_states, n_source, 'uint32');
    dp(1) = 0;

    for src = 1:n_source
        next_dp = -inf(num_states, 1);
        for mask_value = 0:(num_states - 1)
            mask = uint32(mask_value);
            state_idx = mask_value + 1;
            if ~isfinite(dp(state_idx))
                continue;
            end

            % Skip this source frame.
            if dp(state_idx) > next_dp(state_idx)
                next_dp(state_idx) = dp(state_idx);
                parent_col(state_idx, src) = uint16(0);
                parent_state(state_idx, src) = mask;
            end

            % Assign this source frame to one still-unmatched GT cell.
            for col = 1:n_target
                if ~bitget(mask, col)
                    new_mask = bitset(mask, col);
                    new_idx = double(new_mask) + 1;
                    candidate = dp(state_idx) + score_matrix(src, col);
                    if candidate > next_dp(new_idx)
                        next_dp(new_idx) = candidate;
                        parent_col(new_idx, src) = uint16(col);
                        parent_state(new_idx, src) = mask;
                    end
                end
            end
        end
        dp = next_dp;
    end

    assignment = zeros(n_target, 1);
    mask = uint32(num_states - 1);
    if ~isfinite(dp(double(mask) + 1))
        error('Unable to assign all target cells.');
    end

    for src = n_source:-1:1
        idx = double(mask) + 1;
        col = double(parent_col(idx, src));
        prev_mask = parent_state(idx, src);
        if col > 0
            assignment(col) = src;
        end
        mask = prev_mask;
        if mask == 0 && all(assignment > 0)
            break;
        end
    end

    if any(assignment == 0)
        error('Assignment reconstruction failed.');
    end
end

function v = safeMedianPositive(values)
    values = values(values > 0 & isfinite(values));
    if isempty(values)
        v = 1;
    else
        v = median(values);
    end
end

function [aligned_patch, shift_yx] = localRegisterPatch(source_patch, target_patch, target_mask, margin)
    source_patch = double(source_patch);
    target_patch = double(target_patch);
    target_mask = logical(target_mask);

    template = target_patch;
    template(~target_mask) = 0;
    template = normalizeForNcc(template);

    search = padarray(source_patch, [margin, margin], 0, 'both');
    search = normalizeForNcc(search);
    c = normxcorr2(template, search);
    [~, imax] = max(c(:));
    [ypeak, xpeak] = ind2sub(size(c), imax);
    dy = ypeak - size(template, 1) - margin;
    dx = xpeak - size(template, 2) - margin;

    H = size(source_patch, 1);
    W = size(source_patch, 2);
    padded = padarray(source_patch, [margin, margin], 0, 'both');
    rows = (1:H) + margin + dy;
    cols = (1:W) + margin + dx;
    rows = max(1, min(size(padded, 1), rows));
    cols = max(1, min(size(padded, 2), cols));
    aligned_patch = padded(rows, cols);
    shift_yx = [dy, dx];
end

function patch = cleanWithMask(patch, mask, cfg)
    patch = double(patch);
    mask = logical(mask);
    if cfg.keep_only_gt_mask
        mask_out = expandedMask(mask, cfg.output_mask_dilate_px);
        patch(~mask_out) = 0;
    end
    patch = robustZeroPatch(patch, mask, cfg.background_percentile);
end

function mask_out = expandedMask(mask, radius_px)
    mask_out = logical(mask);
    if nargin < 2 || radius_px <= 0 || ~any(mask_out(:))
        return;
    end
    if exist('strel', 'file') == 2 && exist('imdilate', 'file') == 2
        mask_out = imdilate(mask_out, strel('disk', radius_px));
    else
        kernel = true(2 * radius_px + 1);
        mask_out = conv2(double(mask_out), double(kernel), 'same') > 0;
    end
end

function patch = robustZeroPatch(patch, mask, pct)
    patch = double(patch);
    if nargin < 2 || isempty(mask)
        mask = true(size(patch));
    end
    vals = patch(mask & isfinite(patch));
    if isempty(vals)
        patch(:) = 0;
        return;
    end
    offset = prctile(vals, pct);
    patch(mask) = patch(mask) - offset;
    patch(~isfinite(patch)) = 0;
end

function out = normalizeForNcc(in)
    out = double(in);
    out(~isfinite(out)) = 0;
    out = out - mean(out(:));
    sd = std(out(:));
    if sd > 0
        out = out / sd;
    end
end

function out = normalize01(in)
    in = double(in);
    mn = min(in(:));
    mx = max(in(:));
    out = (in - mn) / (mx - mn + eps);
end

function patch = cropCenteredWithPadding(img, cy, cx, H, W, pad_value)
    pad_y = ceil(H / 2) + 2;
    pad_x = ceil(W / 2) + 2;
    img_p = padarray(img, [pad_y, pad_x], pad_value, 'both');
    cy_p = round(cy) + pad_y;
    cx_p = round(cx) + pad_x;
    rows = cy_p - floor((H - 1) / 2) : cy_p + ceil((H - 1) / 2);
    cols = cx_p - floor((W - 1) / 2) : cx_p + ceil((W - 1) / 2);
    patch = img_p(rows, cols);
end

function printAssignmentReport(method_name, report, gt)
    fprintf('     assignment for %s (source -> canonical GT):\n', upper(method_name));
    for row = 1:size(report.assignment, 1)
        src = report.assignment(row, 1);
        dst = report.assignment(row, 2);
        fprintf('       src %02d -> gt %02d (%s), score %.4f\n', ...
            src, dst, gt.label{dst}, report.score_matrix(src, dst));
    end
end

function ovd_report = computeOvdReport(aligned, method_specs, cfg)
    gt = aligned.gt;
    gt_params = gt.physical_params;
    aligned_dx_um = gt_params.dx;
    pixel_area_um2 = aligned_dx_um^2;
    n_cells = size(gt.stack, 3);

    gt_ovd = zeros(n_cells, 1);
    gt_volume = zeros(n_cells, 1);
    for c = 1:n_cells
        mask = expandedMask(gt.mask_stack(:, :, c), cfg.ovd_mask_dilate_px);
        gt_ovd(c) = phaseToOvd(gt.stack(:, :, c), mask, gt_params.lambda, pixel_area_um2);
        gt_volume(c) = gt_ovd(c) / gt_params.delta_n;
    end

    ovd_report = struct([]);
    row = 0;
    for m = 1:numel(method_specs)
        method_name = method_specs(m).name;
        if ~isfield(aligned, method_name)
            continue;
        end

        params = aligned.(method_name).method_params;
        for c = 1:n_cells
            row = row + 1;
            mask = expandedMask(gt.mask_stack(:, :, c), cfg.ovd_mask_dilate_px);
            method_ovd = phaseToOvd(aligned.(method_name).ph_stack(:, :, c), ...
                mask, params.lambda, pixel_area_um2);
            method_volume = method_ovd / params.delta_n;

            ovd_report(row).method = method_name; %#ok<AGROW>
            ovd_report(row).cell_index = c; %#ok<AGROW>
            ovd_report(row).ovd_group = gt.label{c}; %#ok<AGROW>
            ovd_report(row).source_frame = aligned.(method_name).selected_source_frames(c); %#ok<AGROW>
            ovd_report(row).lambda_um = params.lambda; %#ok<AGROW>
            ovd_report(row).native_dx_um = params.dx; %#ok<AGROW>
            ovd_report(row).aligned_dx_um = aligned_dx_um; %#ok<AGROW>
            ovd_report(row).polymerRI = params.polymerRI; %#ok<AGROW>
            ovd_report(row).n_immersion = params.n_immersion; %#ok<AGROW>
            ovd_report(row).delta_n = params.delta_n; %#ok<AGROW>
            ovd_report(row).gt_ovd_um3 = gt_ovd(c); %#ok<AGROW>
            ovd_report(row).method_ovd_um3 = method_ovd; %#ok<AGROW>
            ovd_report(row).ovd_relative_error = safeRelativeError(method_ovd, gt_ovd(c)); %#ok<AGROW>
            ovd_report(row).gt_volume_um3 = gt_volume(c); %#ok<AGROW>
            ovd_report(row).method_volume_um3 = method_volume; %#ok<AGROW>
            ovd_report(row).volume_relative_error = safeRelativeError(method_volume, gt_volume(c)); %#ok<AGROW>
        end
    end
end

function rel = safeRelativeError(value, reference)
    if ~isfinite(reference) || abs(reference) < eps
        rel = NaN;
    else
        rel = (value - reference) / reference;
    end
end

function ovd = phaseToOvd(phase_img, mask, lambda_um, pixel_area_um2)
    values = double(phase_img(mask));
    values = values(isfinite(values));
    ovd = sum(values) * (lambda_um / (2 * pi)) * pixel_area_um2;
end

function writeOvdCsv(ovd_report, out_path)
    if isempty(ovd_report)
        return;
    end
    T = struct2table(ovd_report);
    writetable(T, out_path);
end

function makeOvdPlots(ovd_report, cfg)
    if isempty(ovd_report)
        fprintf('  -> OVD plots skipped: empty report.\n');
        return;
    end

    methods = unique({ovd_report.method}, 'stable');
    cells = unique([ovd_report.cell_index]);
    err_matrix = NaN(numel(methods), numel(cells));
    for i = 1:numel(ovd_report)
        m_idx = find(strcmp(methods, ovd_report(i).method), 1);
        c_idx = find(cells == ovd_report(i).cell_index, 1);
        err_matrix(m_idx, c_idx) = ovd_report(i).ovd_relative_error * 100;
    end

    fig_line = figure('Name', 'OVD Relative Error Lines', 'Color', 'w', 'Position', [100 100 1200 520]);
    hold on; grid on; box on;
    plot([1, numel(cells)], [0, 0], 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    for m = 1:numel(methods)
        plot(1:numel(cells), err_matrix(m, :), '-o', 'LineWidth', 1.6, ...
            'MarkerSize', 5, 'DisplayName', upper(methods{m}));
    end
    cell_labels = arrayfun(@(c) sprintf('C%d', c), cells, 'UniformOutput', false);
    set(gca, 'XTick', 1:numel(cells), 'XTickLabel', cell_labels, 'XTickLabelRotation', 45);
    ylabel('OVD relative error [%]', 'FontWeight', 'bold');
    title('OVD Relative Error by Cell', 'FontWeight', 'bold');
    legend('Location', 'eastoutside', 'Interpreter', 'none');
    saveas(fig_line, fullfile(cfg.output_dir, 'ovd_relative_error_lines.png'));
    close(fig_line);

    fig_box = figure('Name', 'OVD Relative Error Boxplot', 'Color', 'w', 'Position', [150 150 1000 520]);
    valid_vals = [];
    valid_groups = {};
    for m = 1:numel(methods)
        vals = err_matrix(m, :);
        vals = vals(isfinite(vals));
        valid_vals = [valid_vals, vals]; %#ok<AGROW>
        valid_groups = [valid_groups, repmat(methods(m), 1, numel(vals))]; %#ok<AGROW>
    end
    if ~isempty(valid_vals) && exist('boxplot', 'file') == 2
        boxplot(valid_vals, valid_groups, 'LabelOrientation', 'inline');
        grid on; box on;
    else
        hold on; grid on; box on;
        for m = 1:numel(methods)
            scatter(m * ones(1, size(err_matrix, 2)), err_matrix(m, :), 30, 'filled');
        end
        set(gca, 'XTick', 1:numel(methods), 'XTickLabel', upper(methods), 'XTickLabelRotation', 45);
    end
    ylabel('OVD relative error [%]', 'FontWeight', 'bold');
    title('OVD Relative Error Distribution', 'FontWeight', 'bold');
    saveas(fig_box, fullfile(cfg.output_dir, 'ovd_relative_error_boxplot.png'));
    close(fig_box);

    mean_err = NaN(1, numel(methods));
    for m = 1:numel(methods)
        vals = err_matrix(m, :);
        mean_err(m) = mean(vals(isfinite(vals)));
    end
    fig_bar = figure('Name', 'Mean OVD Relative Error', 'Color', 'w', 'Position', [200 200 1000 520]);
    bar(mean_err, 'FaceColor', [0.2 0.45 0.8], 'EdgeColor', 'k');
    grid on; box on;
    hold on;
    plot(xlim, [0, 0], 'k-', 'LineWidth', 1.2);
    set(gca, 'XTick', 1:numel(methods), 'XTickLabel', upper(methods), 'XTickLabelRotation', 45);
    ylabel('Mean OVD relative error [%]', 'FontWeight', 'bold');
    title('Mean OVD Relative Error by Method', 'FontWeight', 'bold');
    for m = 1:numel(methods)
        if isfinite(mean_err(m))
            text(m, mean_err(m), sprintf(' %.2f%%', mean_err(m)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', valueLabelAlignment(mean_err(m)), ...
                'FontWeight', 'bold');
        end
    end
    saveas(fig_bar, fullfile(cfg.output_dir, 'ovd_mean_relative_error_bar.png'));
    close(fig_bar);
end

function align = valueLabelAlignment(value)
    if value >= 0
        align = 'bottom';
    else
        align = 'top';
    end
end

function makeVerificationGallery(aligned, method_specs, cfg)
    method_names = {};
    for i = 1:numel(method_specs)
        if isfield(aligned, method_specs(i).name)
            method_names{end + 1} = method_specs(i).name; %#ok<AGROW>
        end
    end
    if isempty(method_names)
        fprintf('  -> No method stacks available; gallery not generated.\n');
        return;
    end

    for m = 1:numel(method_names)
        name = method_names{m};
        stack = aligned.(name).ph_stack;
        fig = figure('Name', sprintf('%s aligned verification', upper(name)), ...
            'Position', [100 100 1800 900], 'Color', 'w');
        tlo = tiledlayout(4, 9, 'TileSpacing', 'none', 'Padding', 'compact');
        title(tlo, sprintf('%s aligned to canonical GT: rows 1/3 method, rows 2/4 GT', upper(name)), ...
            'FontSize', 15, 'FontWeight', 'bold');

        for i = 1:9
            nexttile(i);
            imagesc(stack(:, :, i)); axis image off; colormap jet;
            title(sprintf('6k %s C%d', upper(name), i), 'FontSize', 8);
        end
        for i = 1:9
            nexttile(i + 9);
            imagesc(aligned.gt.stack(:, :, i)); axis image off; colormap jet;
            title(sprintf('GT C%d', i), 'FontSize', 8);
        end
        for i = 10:18
            nexttile(i + 9);
            imagesc(stack(:, :, i)); axis image off; colormap jet;
            title(sprintf('15k %s C%d', upper(name), i - 9), 'FontSize', 8);
        end
        for i = 10:18
            nexttile(i + 18);
            imagesc(aligned.gt.stack(:, :, i)); axis image off; colormap jet;
            title(sprintf('GT C%d', i - 9), 'FontSize', 8);
        end

        out_png = fullfile(cfg.output_dir, sprintf('%s_alignment_verification_gallery.png', name));
        saveas(fig, out_png);
        close(fig);
        fprintf('  -> Gallery saved: %s\n', out_png);
    end
end
