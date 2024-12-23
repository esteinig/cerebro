// Provides access to system resources and compiler functions

use std::{
    borrow::Borrow,
    collections::HashMap,
    path::{Path, PathBuf},
    sync::OnceLock,
};

use chrono::{DateTime, Datelike, Local};
use parking_lot::Mutex;
use send_wrapper::SendWrapper;
use typst::{
    diag::{EcoString, FileError, FileResult},
    foundations::{Bytes, Datetime},
    model::Document,
    syntax::{FileId, Source, VirtualPath},
    text::{Font, FontBook},
    utils::LazyHash,
    Library, World,
};

use wasm_bindgen::JsValue;

use std::cell::OnceCell;


pub struct FileEntry {
    bytes: OnceCell<Bytes>,
    source: Source,
}

impl FileEntry {
    pub fn new(id: FileId, text: String) -> Self {
        Self {
            bytes: OnceCell::new(),
            source: Source::new(id, text),
        }
    }

    pub fn source(&self) -> Source {
        self.source.clone()
    }

    pub fn bytes(&self) -> Bytes {
        self.bytes
            .get_or_init(|| Bytes::from(self.source.text().as_bytes()))
            .clone()
    }
}

pub struct EmbeddedFonts {
    pub regular: &'static [u8],
    pub bold: &'static [u8],
    pub italic: &'static [u8],
}

pub fn get_embedded_fonts() -> EmbeddedFonts {
    EmbeddedFonts {
        regular: include_bytes!("../fonts/IBMPlexSans-Regular.ttf"),
        bold: include_bytes!("../fonts/IBMPlexSans-Bold.ttf"),
        italic: include_bytes!("../fonts/IBMPlexSans-Italic.ttf"),
    }
}

pub struct SystemWorld {
    /// The root relative to which absolute paths are resolved.
    root: PathBuf,
    /// The input source.
    main: FileId,
    /// Typst's standard library.
    library: LazyHash<Library>,
    /// Metadata about discovered fonts.
    book: LazyHash<FontBook>,
    /// Storage of fonts
    fonts: Vec<Font>,

    files: Mutex<HashMap<FileId, FileEntry>>,
    now: OnceLock<DateTime<Local>>,
    request_data: SendWrapper<js_sys::Function>,

}

impl SystemWorld {
    pub fn new(root: String, request_data: &js_sys::Function) -> SystemWorld {

        let (book, fonts) = SystemWorld::start_embedded_fonts();

        Self {
            root: PathBuf::from(root),
            main: FileId::new(None, VirtualPath::new("")),
            library: LazyHash::new(Library::default()),
            book: LazyHash::new(book),
            fonts,
            files: Mutex::default(),
            now: OnceLock::new(),
            request_data: SendWrapper::new(request_data.clone()),
        }
    }

    pub fn compile(&mut self, text: String, path: String) -> Result<Document, JsValue> {

        self.reset();

        self.main = FileId::new(None, VirtualPath::new(path));

        self.files
            .get_mut()
            .insert(self.main, FileEntry::new(self.main, text));

        typst::compile(self)
            .output
            .map_err(|_| JsValue::from("Error in compiler"))
    }

    pub fn add_font(&mut self, data: Vec<u8>) {

        let buffer = Bytes::from(data);
        let mut font_infos = Vec::new();
        for font in Font::iter(buffer) {
            font_infos.push(font.info().clone());
            self.fonts.push(font)
        }
        if font_infos.len() > 0 {
            for info in font_infos {
                self.book.push(info);
            }
        }
    }

    fn reset(&mut self) {
        self.files.get_mut().clear();
        self.now.take();
    }

    fn request_data(&self, param1: String) -> Result<JsValue, JsValue> {
        return self.request_data.call1(&JsValue::NULL, &param1.into());
    }

    fn read_file(&self, path: &Path) -> FileResult<String> {
        
        let f = |e: JsValue| {
            if let Some(value) = e.as_f64() {
                return match value as i64 {
                    2 => FileError::NotFound(path.to_path_buf()),
                    3 => FileError::AccessDenied,
                    4 => FileError::IsDirectory,
                    _ => FileError::Other(Some(EcoString::from("see console for details"))),
                };
            }
            FileError::Other(e.as_string().map(EcoString::from))
        };

        Ok(self
            .request_data(path.to_str().unwrap().to_owned())
            .map_err(f)?
            .as_string()
            .unwrap())
    }


    fn file_entry<F, T>(&self, id: FileId, f: F) -> FileResult<T>
    where
        F: FnOnce(&mut FileEntry) -> T,
    {
        let mut map = self.files.lock();
        if !map.contains_key(&id) {
            let path = self.root.clone();

            let text =
                self.read_file(&id.vpath().resolve(&path).ok_or(FileError::AccessDenied)?)?;

            map.insert(id, FileEntry::new(id, text));
        }
        Ok(f(map.get_mut(&id).unwrap()))
    }

    fn start_embedded_fonts() -> (FontBook, Vec<Font>) {
        let mut book = FontBook::new();
        let mut fonts = Vec::new();

        for data in typst_assets::fonts() {
            let buffer = Bytes::from_static(data);
            for font in Font::iter(buffer) {
                book.push(font.info().clone());
                fonts.push(font);
            }
        }
        let embedded_fonts = get_embedded_fonts();
        
        for data in [embedded_fonts.regular, embedded_fonts.bold, embedded_fonts.italic] {
            let buffer = Bytes::from_static(data);
            for font in Font::iter(buffer) {
                book.push(font.info().clone());
                fonts.push(font);
            }
        }

        return (book, fonts);
    }
}

impl World for SystemWorld {
    fn library(&self) -> &LazyHash<Library> {
        &self.library
    }

    fn book(&self) -> &LazyHash<FontBook> {
        &self.book
    }

    fn main(&self) -> FileId {
        self.main
    }

    fn source(&self, id: FileId) -> FileResult<Source> {
        self.file_entry(id, |f| f.source())
    }

    fn file(&self, id: FileId) -> FileResult<Bytes> {
        self.file_entry(id, |f| f.bytes())
    }

    fn font(&self, index: usize) -> Option<Font> {
        Some(self.fonts[index].clone())
    }

    fn today(&self, offset: Option<i64>) -> Option<Datetime> {
        let now = self.now.get_or_init(chrono::Local::now);

        let naive = match offset {
            None => now.naive_local(),
            Some(o) => now.naive_utc() + chrono::Duration::hours(o),
        };

        Datetime::from_ymd(
            naive.year(),
            naive.month().try_into().ok()?,
            naive.day().try_into().ok()?,
        )
    }
}